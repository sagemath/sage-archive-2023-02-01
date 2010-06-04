r"""
A Collection of Constructors of Common Graphs

Usage
=====

To see a list of all graph constructors, type "graphs." and then
press the tab key. The documentation for each constructor includes
information about each graph, which provides a useful reference.


Plotting
========

All graphs (i.e., networks) have an associated Sage
graphics object, which you can display::

    sage: G = graphs.WheelGraph(15)
    sage: P = G.plot()
    sage: P.show() # long time

If you create a graph in Sage using the ``Graph``
command, then plot that graph, the positioning of nodes is
determined using the spring-layout algorithm. For the special graph
constructors, which you get using ``graphs.[tab]``, the
positions are preset. For example, consider the Petersen graph with
default node positioning vs. the Petersen graph constructed by this
database::

    sage: petersen_spring = Graph({0:[1,4,5], 1:[0,2,6], 2:[1,3,7], 3:[2,4,8], 4:[0,3,9], 5:[0,7,8], 6:[1,8,9], 7:[2,5,9], 8:[3,5,6], 9:[4,6,7]})
    sage: petersen_spring.show() # long time
    sage: petersen_database = graphs.PetersenGraph()
    sage: petersen_database.show() # long time

For all the constructors in this database (except the octahedral,
dodecahedral, random and empty graphs), the position dictionary is
filled in, instead of using the spring-layout algorithm.

For further visual examples and explanation, see the docstrings
below, particularly for
:meth:`CycleGraph <GraphGenerators.CycleGraph>`,
:meth:`StarGraph <GraphGenerators.StarGraph>`,
:meth:`WheelGraph <GraphGenerators.WheelGraph>`,
:meth:`CompleteGraph <GraphGenerators.CompleteGraph>`, and
:meth:`CompleteBipartiteGraph <GraphGenerators.CompleteBipartiteGraph>`.


.. _organization:

Organization
============

The constructors available in this database are
organized as follows.


Basic structures
----------------

- :meth:`BarbellGraph <GraphGenerators.BarbellGraph>`
- :meth:`BullGraph <GraphGenerators.BullGraph>`
- :meth:`CircularLadderGraph <GraphGenerators.CircularLadderGraph>`
- :meth:`ClawGraph <GraphGenerators.ClawGraph>`
- :meth:`CycleGraph <GraphGenerators.CycleGraph>`
- :meth:`DiamondGraph <GraphGenerators.DiamondGraph>`
- :meth:`EmptyGraph <GraphGenerators.EmptyGraph>`
- :meth:`Grid2dGraph <GraphGenerators.Grid2dGraph>`
- :meth:`GridGraph <GraphGenerators.GridGraph>`
- :meth:`HouseGraph <GraphGenerators.HouseGraph>`
- :meth:`HouseXGraph <GraphGenerators.HouseXGraph>`
- :meth:`KrackhardtKiteGraph <GraphGenerators.KrackhardtKiteGraph>`
- :meth:`LadderGraph <GraphGenerators.LadderGraph>`
- :meth:`LCFGraph <GraphGenerators.LCFGraph>`
- :meth:`LollipopGraph <GraphGenerators.LollipopGraph>`
- :meth:`PathGraph <GraphGenerators.PathGraph>`
- :meth:`StarGraph <GraphGenerators.StarGraph>`
- :meth:`ToroidalGrid2dGraph <GraphGenerators.ToroidalGrid2dGraph>`
- :meth:`WheelGraph <GraphGenerators.WheelGraph>`


Platonic solids
---------------

- :meth:`DodecahedralGraph <GraphGenerators.DodecahedralGraph>`
- :meth:`HexahedralGraph <GraphGenerators.HexahedralGraph>`
- :meth:`IcosahedralGraph <GraphGenerators.IcosahedralGraph>`
- :meth:`OctahedralGraph <GraphGenerators.OctahedralGraph>`
- :meth:`TetrahedralGraph <GraphGenerators.TetrahedralGraph>`


Named Graphs
------------

- :meth:`ChvatalGraph <GraphGenerators.ChvatalGraph>`
- :meth:`DesarguesGraph <GraphGenerators.DesarguesGraph>`
- :meth:`FlowerSnark <GraphGenerators.FlowerSnark>`
- :meth:`FruchtGraph <GraphGenerators.FruchtGraph>`
- :meth:`HeawoodGraph <GraphGenerators.HeawoodGraph>`
- :meth:`HigmanSimsGraph <GraphGenerators.HigmanSimsGraph>`
- :meth:`HoffmanSingletonGraph <GraphGenerators.HoffmanSingletonGraph>`
- :meth:`MoebiusKantorGraph <GraphGenerators.MoebiusKantorGraph>`
- :meth:`PappusGraph <GraphGenerators.PappusGraph>`
- :meth:`PetersenGraph <GraphGenerators.PetersenGraph>`
- :meth:`ThomsenGraph <GraphGenerators.ThomsenGraph>`


Families of graphs
------------------

- :meth:`BalancedTree <GraphGenerators.BalancedTree>`
- :meth:`BubbleSortGraph <GraphGenerators.BubbleSortGraph>`
- :meth:`CirculantGraph <GraphGenerators.CirculantGraph>`
- :meth:`CompleteBipartiteGraph <GraphGenerators.CompleteBipartiteGraph>`
- :meth:`CompleteGraph <GraphGenerators.CompleteGraph>`
- :meth:`CubeGraph <GraphGenerators.CubeGraph>`
- :meth:`FibonacciTree <GraphGenerators.FibonacciTree>`
- :meth:`FuzzyBallGraph <GraphGenerators.FuzzyBallGraph>`
- :meth:`GeneralizedPetersenGraph <GraphGenerators.GeneralizedPetersenGraph>`
- :meth:`HanoiTowerGraph <GraphGenerators.HanoiTowerGraph>`
- :meth:`HyperStarGraph <GraphGenerators.HyperStarGraph>`
- :meth:`KneserGraph <GraphGenerators.KneserGraph>`
- :meth:`LCFGraph <GraphGenerators.LCFGraph>`
- :meth:`NKStarGraph <GraphGenerators.NKStarGraph>`
- :meth:`NStarGraph <GraphGenerators.NStarGraph>`
- :meth:`OddGraph <GraphGenerators.OddGraph>`
- :meth:`trees <GraphGenerators.trees>`


Pseudofractal graphs
--------------------

- :meth:`DorogovtsevGoltsevMendesGraph <GraphGenerators.DorogovtsevGoltsevMendesGraph>`


Random graphs
-------------

- :meth:`RandomBarabasiAlbert <GraphGenerators.RandomBarabasiAlbert>`
- :meth:`RandomBipartite <GraphGenerators.RandomBipartite>`
- :meth:`RandomGNM <GraphGenerators.RandomGNM>`
- :meth:`RandomGNP <GraphGenerators.RandomGNP>`
- :meth:`RandomHolmeKim <GraphGenerators.RandomHolmeKim>`
- :meth:`RandomInterval <GraphGenerators.RandomInterval>`
- :meth:`RandomLobster <GraphGenerators.RandomLobster>`
- :meth:`RandomNewmanWattsStrogatz <GraphGenerators.RandomNewmanWattsStrogatz>`
- :meth:`RandomRegular <GraphGenerators.RandomRegular>`
- :meth:`RandomShell <GraphGenerators.RandomShell>`
- :meth:`RandomTreePowerlaw <GraphGenerators.RandomTreePowerlaw>`


Graphs with a given degree sequence
-----------------------------------

- :meth:`DegreeSequence <GraphGenerators.DegreeSequence>`
- :meth:`DegreeSequenceBipartite <GraphGenerators.DegreeSequenceBipartite>`
- :meth:`DegreeSequenceConfigurationModel <GraphGenerators.DegreeSequenceConfigurationModel>`
- :meth:`DegreeSequenceExpected <GraphGenerators.DegreeSequenceExpected>`
- :meth:`DegreeSequenceTree <GraphGenerators.DegreeSequenceTree>`


Miscellaneous
-------------

- :meth:`WorldMap <GraphGenerators.WorldMap>`


AUTHORS:

- Robert Miller (2006-11-05): initial version, empty, random, petersen

- Emily Kirkman (2006-11-12): basic structures, node positioning for
  all constructors

- Emily Kirkman (2006-11-19): docstrings, examples

- William Stein (2006-12-05): Editing.

- Robert Miller (2007-01-16): Cube generation and plotting

- Emily Kirkman (2007-01-16): more basic structures, docstrings

- Emily Kirkman (2007-02-14): added more named graphs

- Robert Miller (2007-06-08-11): Platonic solids, random graphs,
  graphs with a given degree sequence, random directed graphs

- Robert Miller (2007-10-24): Isomorph free exhaustive generation

- Nathann Cohen (2009-08-12): WorldMap

- Michael Yurko (2009-9-01): added hyperstar, (n,k)-star, n-star, and
  bubblesort graphs

- Anders Jonsson (2009-10-15): added generalized Petersen graphs

- Harald Schilly and Yann Laigle-Chapuy (2010-03-24): added Fibonacci Tree

- Jason Grout (2010-06-04): cospectral_graphs
"""

###########################################################################

#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#           Copyright (C) 2009 Michael C. Yurko <myurko@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
###########################################################################

# import from Python standard library
from math import sin, cos, pi

# import from Sage library
import graph
from sage.misc.randstate import current_randstate

class GraphGenerators():
    r"""
    A class consisting of constructors for several common graphs, as
    well as orderly generation of isomorphism class representatives. See the
    section :ref:`organization` for a list of supported constructors.

    A list of all graphs and graph structures (other than isomorphism class
    representatives) in this database is available via tab completion. Type
    "graphs." and then hit the tab key to see which graphs are available.

    The docstrings include educational information about each named
    graph with the hopes that this class can be used as a reference.

    For all the constructors in this class (except the octahedral,
    dodecahedral, random and empty graphs), the position dictionary is
    filled to override the spring-layout algorithm.


    ORDERLY GENERATION::

        graphs(vertices, property=lambda x: True, augment='edges', size=None)

    This syntax accesses the generator of isomorphism class
    representatives. Iterates over distinct, exhaustive
    representatives.

    INPUT:

    - ``vertices`` -- natural number.

    - ``property`` -- (default: ``lambda x: True``) any property to be tested
      on graphs before generation. (Ignored if ``deg_seq`` is specified.)

    - ``augment`` -- (default: ``'edges'``) possible values:

      - ``'vertices'`` -- augments by adding a vertex and
        edges incident to that vertex. In this case, all graphs up to
        ``n=vertices`` are generated. If for any graph G satisfying the
        property, every subgraph, obtained from G by deleting one vertex
        and only edges incident to that vertex, satisfies the property,
        then this will generate all graphs with that property. If this does
        not hold, then all the graphs generated will satisfy the property,
        but there will be some missing.

      - ``'edges'`` -- augments a fixed number of vertices by
        adding one edge. In this case, all graphs on exactly ``n=vertices`` are
        generated. If for any graph G satisfying the property, every
        subgraph, obtained from G by deleting one edge but not the vertices
        incident to that edge, satisfies the property, then this will
        generate all graphs with that property. If this does not hold, then
        all the graphs generated will satisfy the property, but there will
        be some missing.

    - ``size`` -- (default: ``None``) the size of the graph to be generated.

    - ``deg_seq`` -- (default: ``None``) a sequence of non-negative integers,
      or ``None``. If specified, the generated graphs will have these
      integers for degrees. In this case, property and size are both
      ignored.

    - ``loops`` -- (default: ``False``) whether to allow loops in the graph
      or not.

    - ``implementation`` -- (default: ``'c_graph'``) which underlying
      implementation to use (see ``Graph?``).

    - ``sparse`` -- (default: ``True``) ignored if implementation is not
      ``'c_graph'``.

    EXAMPLES:

    Print graphs on 3 or less vertices::

        sage: for G in graphs(3, augment='vertices'):
        ...    print G
        Graph on 0 vertices
        Graph on 1 vertex
        Graph on 2 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 2 vertices
        Graph on 3 vertices

    Note that we can also get graphs with underlying Cython implementation::

        sage: for G in graphs(3, augment='vertices', implementation='c_graph'):
        ...    print G
        Graph on 0 vertices
        Graph on 1 vertex
        Graph on 2 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 2 vertices
        Graph on 3 vertices

    Print graphs on 3 vertices.

    ::

        sage: for G in graphs(3):
        ...    print G
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices

    Generate all graphs with 5 vertices and 4 edges.

    ::

        sage: L = graphs(5, size=4)
        sage: len(list(L))
        6

    Generate all graphs with 5 vertices and up to 4 edges.

    ::

        sage: L = list(graphs(5, lambda G: G.size() <= 4))
        sage: len(L)
        14
        sage: graphs_list.show_graphs(L) # long time

    Generate all graphs with up to 5 vertices and up to 4 edges.

    ::

        sage: L = list(graphs(5, lambda G: G.size() <= 4, augment='vertices'))
        sage: len(L)
        31
        sage: graphs_list.show_graphs(L)              # long time

    Generate all graphs with degree at most 2, up to 6 vertices.

    ::

        sage: property = lambda G: ( max([G.degree(v) for v in G] + [0]) <= 2 )
        sage: L = list(graphs(6, property, augment='vertices'))
        sage: len(L)
        45

    Generate all bipartite graphs on up to 7 vertices: (see
    http://www.research.att.com/~njas/sequences/A033995)

    ::

        sage: L = list( graphs(7, lambda G: G.is_bipartite(), augment='vertices') )
        sage: [len([g for g in L if g.order() == i]) for i in [1..7]]
        [1, 2, 3, 7, 13, 35, 88]

    Generate all bipartite graphs on exactly 7 vertices::

        sage: L = list( graphs(7, lambda G: G.is_bipartite()) )
        sage: len(L)
        88

    Generate all bipartite graphs on exactly 8 vertices::

        sage: L = list( graphs(8, lambda G: G.is_bipartite()) ) # long time
        sage: len(L)                                            # long time
        303

    Generate graphs on the fly: (see
    http://www.research.att.com/~njas/sequences/A000088)

    ::

        sage: for i in range(0, 7):
        ...    print len(list(graphs(i)))
        1
        1
        2
        4
        11
        34
        156

    Generate all simple graphs, allowing loops: (see
    http://www.research.att.com/~njas/sequences/A000666)

    ::

        sage: L = list(graphs(5,augment='vertices',loops=True))               # long time
        sage: for i in [0..5]: print i, len([g for g in L if g.order() == i]) # long time
        0 1
        1 2
        2 6
        3 20
        4 90
        5 544

    Generate all graphs with a specified degree sequence (see
    http://www.research.att.com/~njas/sequences/A002851)::

        sage: for i in [4,6,8]:
        ...       print i, len([g for g in graphs(i, deg_seq=[3]*i) if g.is_connected()])
        4 1
        6 2
        8 5
        sage: for i in [4,6,8]:                                                                          # long time
        ...       print i, len([g for g in graphs(i,augment='vertices',deg_seq=[3]*i) if g.is_connected()]) # long time
        4 1
        6 2
        8 5

    ::

        sage: print 10, len([g for g in graphs(10,deg_seq=[3]*10) if g.is_connected()]) # not tested
        10 19

    REFERENCE:

    - Brendan D. McKay, Isomorph-Free Exhaustive generation.  *Journal
      of Algorithms*, Volume 26, Issue 2, February 1998, pages 306-324.
    """

    #######################################################################
    #   Basic Structures
    #######################################################################

    def BarbellGraph(self, n1, n2):
        r"""
        Returns a barbell graph with ``2*n1 + n2`` nodes. The argument ``n1``
        must be greater than or equal to 2.

        A barbell graph is a basic structure that consists of a path graph
        of order ``n2`` connecting two complete graphs of order ``n1`` each.

        This constructor depends on `NetworkX <http://networkx.lanl.gov>`_
        numeric labels. In this case, the ``n1``-th node connects to the
        path graph from one complete graph and the ``n1 + n2 + 1``-th node
        connects to the path graph from the other complete graph.

        INPUT:

        - ``n1`` -- integer `\geq 2`. The order of each of the two
          complete graphs.

        - ``n2`` -- nonnegative integer. The order of the path graph
          connecting the two complete graphs.

        OUTPUT:

        A barbell graph of order ``2*n1 + n2``. A ``ValueError`` is
        returned if ``n1 < 2`` or ``n2 < 0``.

        ALGORITHM:

        Uses `NetworkX <http://networkx.lanl.gov>`_.

        PLOTTING:

        Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, each barbell
        graph will be displayed with the two complete graphs in the
        lower-left and upper-right corners, with the path graph connecting
        diagonally between the two. Thus the ``n1``-th node will be drawn at a
        45 degree angle from the horizontal right center of the first
        complete graph, and the ``n1 + n2 + 1``-th node will be drawn 45
        degrees below the left horizontal center of the second complete graph.

        EXAMPLES:

        Construct and show a barbell graph ``Bar = 4``, ``Bells = 9``::

            sage: g = graphs.BarbellGraph(9, 4); g
            Barbell graph: Graph on 22 vertices
            sage: g.show() # long time

        An ``n1 >= 2``, ``n2 >= 0`` barbell graph has order ``2*n1 + n2``. It
        has the complete graph on ``n1`` vertices as a subgraph. It also has
        the path graph on ``n2`` vertices as a subgraph. ::

            sage: n1 = randint(2, 2*10^2)
            sage: n2 = randint(0, 2*10^2)
            sage: g = graphs.BarbellGraph(n1, n2)
            sage: v = 2*n1 + n2
            sage: g.order() == v
            True
            sage: K_n1 = graphs.CompleteGraph(n1)
            sage: P_n2 = graphs.PathGraph(n2)
            sage: s_K = g.subgraph_search(K_n1, induced=True)
            sage: s_P = g.subgraph_search(P_n2, induced=True)
            sage: K_n1.is_isomorphic(s_K)
            True
            sage: P_n2.is_isomorphic(s_P)
            True

        Create several barbell graphs in a Sage graphics array::

            sage: g = []
            sage: j = []
            sage: for i in range(6):
            ...       k = graphs.BarbellGraph(i + 2, 4)
            ...       g.append(k)
            ...
            sage: for i in range(2):
            ...       n = []
            ...       for m in range(3):
            ...           n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...       j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time

        TESTS:

        The input ``n1`` must be `\geq 2`::

            sage: graphs.BarbellGraph(1, randint(0, 10^6))
            Traceback (most recent call last):
            ...
            ValueError: Invalid graph description, n1 should be >= 2
            sage: graphs.BarbellGraph(randint(-10^6, 1), randint(0, 10^6))
            Traceback (most recent call last):
            ...
            ValueError: Invalid graph description, n1 should be >= 2

        The input ``n2`` must be `\geq 0`::

            sage: graphs.BarbellGraph(randint(2, 10^6), -1)
            Traceback (most recent call last):
            ...
            ValueError: Invalid graph description, n2 should be >= 0
            sage: graphs.BarbellGraph(randint(2, 10^6), randint(-10^6, -1))
            Traceback (most recent call last):
            ...
            ValueError: Invalid graph description, n2 should be >= 0
            sage: graphs.BarbellGraph(randint(-10^6, 1), randint(-10^6, -1))
            Traceback (most recent call last):
            ...
            ValueError: Invalid graph description, n1 should be >= 2
        """
        # sanity checks
        if n1 < 2:
            raise ValueError("Invalid graph description, n1 should be >= 2")
        if n2 < 0:
            raise ValueError("Invalid graph description, n2 should be >= 0")

        pos_dict = {}

        for i in range(n1):
            x = float(cos((pi / 4) - ((2 * pi) / n1) * i) - (n2 / 2) - 1)
            y = float(sin((pi / 4) - ((2 * pi) / n1) * i) - (n2 / 2) - 1)
            j = n1 - 1 - i
            pos_dict[j] = (x, y)
        for i in range(n1, n1 + n2):
            x = float(i - n1 - (n2 / 2) + 1)
            y = float(i - n1 - (n2 / 2) + 1)
            pos_dict[i] = (x, y)
        for i in range(n1 + n2, (2 * n1) + n2):
            x = float(
                cos((5 * (pi / 4)) + ((2 * pi) / n1) * (i - n1 - n2))
                + (n2 / 2) + 2)
            y = float(
                sin((5 * (pi / 4)) + ((2 * pi) / n1) * (i - n1 - n2))
                + (n2 / 2) + 2)
            pos_dict[i] = (x, y)

        import networkx
        G = networkx.barbell_graph(n1, n2)
        return graph.Graph(G, pos=pos_dict, name="Barbell graph")

    def BullGraph(self):
        r"""
        Returns a bull graph with 5 nodes.

        A bull graph is named for its shape. It's a triangle with horns.
        This constructor depends on `NetworkX <http://networkx.lanl.gov>`_
        numeric labeling. For more information, see this
        `Wikipedia article on the bull graph <http://en.wikipedia.org/wiki/Bull_graph>`_.

        PLOTTING:

        Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the bull graph
        is drawn as a triangle with the first node (0) on the bottom. The
        second and third nodes (1 and 2) complete the triangle. Node 3 is
        the horn connected to 1 and node 4 is the horn connected to node
        2.

        ALGORITHM:

        Uses `NetworkX <http://networkx.lanl.gov>`_.

        EXAMPLES:

        Construct and show a bull graph::

            sage: g = graphs.BullGraph(); g
            Bull graph: Graph on 5 vertices
            sage: g.show() # long time

        The bull graph has 5 vertices and 5 edges. Its radius is 2, its
        diameter 3, and its girth 3. The bull graph is planar with chromatic
        number 3 and chromatic index also 3. ::

            sage: g.order(); g.size()
            5
            5
            sage: g.radius(); g.diameter(); g.girth()
            2
            3
            3
            sage: g.chromatic_number()
            3

        The bull graph has chromatic polynomial `x(x - 2)(x - 1)^3` and
        Tutte polynomial `x^4 + x^3 + x^2 y`. Its characteristic polynomial
        is `x(x^2 - x - 3)(x^2 + x - 1)`, which follows from the definition of
        characteristic polynomials for graphs, i.e. `\det(xI - A)`, where
        `x` is a variable, `A` the adjacency matrix of the graph, and `I`
        the identity matrix of the same dimensions as `A`. ::

            sage: chrompoly = g.chromatic_polynomial()
            sage: bool(expand(x * (x - 2) * (x - 1)^3) == chrompoly)
            True
            sage: charpoly = g.characteristic_polynomial()
            sage: M = g.adjacency_matrix(); M
            [0 1 1 0 0]
            [1 0 1 1 0]
            [1 1 0 0 1]
            [0 1 0 0 0]
            [0 0 1 0 0]
            sage: Id = identity_matrix(ZZ, M.nrows())
            sage: D = x*Id - M
            sage: bool(D.determinant() == charpoly)
            True
            sage: bool(expand(x * (x^2 - x - 3) * (x^2 + x - 1)) == charpoly)
            True
        """
        pos_dict = {0:(0,0), 1:(-1,1), 2:(1,1), 3:(-2,2), 4:(2,2)}
        import networkx
        G = networkx.bull_graph()
        return graph.Graph(G, pos=pos_dict, name="Bull graph")

    def CircularLadderGraph(self, n):
        """
        Returns a circular ladder graph with 2\*n nodes.

        A Circular ladder graph is a ladder graph that is connected at the
        ends, i.e.: a ladder bent around so that top meets bottom. Thus it
        can be described as two parallel cycle graphs connected at each
        corresponding node pair.

        This constructor depends on NetworkX numeric labels.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the circular
        ladder graph is displayed as an inner and outer cycle pair, with
        the first n nodes drawn on the inner circle. The first (0) node is
        drawn at the top of the inner-circle, moving clockwise after that.
        The outer circle is drawn with the (n+1)th node at the top, then
        counterclockwise as well.

        EXAMPLES: Construct and show a circular ladder graph with 26 nodes

        ::

            sage: g = graphs.CircularLadderGraph(13)
            sage: g.show() # long time

        Create several circular ladder graphs in a Sage graphics array

        ::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.CircularLadderGraph(i+3)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time
        """
        pos_dict = {}
        for i in range(n):
            x = float(cos((pi/2) + ((2*pi)/n)*i))
            y = float(sin((pi/2) + ((2*pi)/n)*i))
            pos_dict[i] = [x,y]
        for i in range(n,2*n):
            x = float(2*(cos((pi/2) + ((2*pi)/n)*(i-n))))
            y = float(2*(sin((pi/2) + ((2*pi)/n)*(i-n))))
            pos_dict[i] = (x,y)
        import networkx
        G = networkx.circular_ladder_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Circular Ladder graph")

    def ClawGraph(self):
        """
        Returns a claw graph.

        A claw graph is named for its shape. It is actually a complete
        bipartite graph with (n1, n2) = (1, 3).

        PLOTTING: See CompleteBipartiteGraph.

        EXAMPLES: Show a Claw graph

        ::

            sage: (graphs.ClawGraph()).show() # long time

        Inspect a Claw graph

        ::

            sage: G = graphs.ClawGraph()
            sage: G
            Claw graph: Graph on 4 vertices
        """
        pos_dict = {0:(0,1),1:(-1,0),2:(0,0),3:(1,0)}
        import networkx
        G = networkx.complete_bipartite_graph(1,3)
        return graph.Graph(G, pos=pos_dict, name="Claw graph")

    def CycleGraph(self, n):
        r"""
        Returns a cycle graph with n nodes.

        A cycle graph is a basic structure which is also typically called
        an n-gon.

        This constructor is dependent on vertices numbered 0 through n-1 in
        NetworkX ``cycle_graph()``

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, each cycle
        graph will be displayed with the first (0) node at the top, with
        the rest following in a counterclockwise manner.

        The cycle graph is a good opportunity to compare efficiency of
        filling a position dictionary vs. using the spring-layout algorithm
        for plotting. Because the cycle graph is very symmetric, the
        resulting plots should be similar (in cases of small n).

        Filling the position dictionary in advance adds O(n) to the
        constructor.

        EXAMPLES: Compare plotting using the predefined layout and
        networkx::

            sage: import networkx
            sage: n = networkx.cycle_graph(23)
            sage: spring23 = Graph(n)
            sage: posdict23 = graphs.CycleGraph(23)
            sage: spring23.show() # long time
            sage: posdict23.show() # long time

        We next view many cycle graphs as a Sage graphics array. First we
        use the ``CycleGraph`` constructor, which fills in the
        position dictionary::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.CycleGraph(i+3)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time

        Compare to plotting with the spring-layout algorithm::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    spr = networkx.cycle_graph(i+3)
            ...    k = Graph(spr)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time
        """
        pos_dict = {}
        for i in range(n):
            x = float(cos((pi/2) + ((2*pi)/n)*i))
            y = float(sin((pi/2) + ((2*pi)/n)*i))
            pos_dict[i] = (x,y)
        import networkx
        G = networkx.cycle_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Cycle graph")

    def DiamondGraph(self):
        """
        Returns a diamond graph with 4 nodes.

        A diamond graph is a square with one pair of diagonal nodes
        connected.

        This constructor depends on NetworkX numeric labeling.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the diamond
        graph is drawn as a diamond, with the first node on top, second on
        the left, third on the right, and fourth on the bottom; with the
        second and third node connected.

        EXAMPLES: Construct and show a diamond graph

        ::

            sage: g = graphs.DiamondGraph()
            sage: g.show() # long time
        """
        pos_dict = {0:(0,1),1:(-1,0),2:(1,0),3:(0,-1)}
        import networkx
        G = networkx.diamond_graph()
        return graph.Graph(G, pos=pos_dict, name="Diamond Graph")

    def EmptyGraph(self):
        """
        Returns an empty graph (0 nodes and 0 edges).

        This is useful for constructing graphs by adding edges and vertices
        individually or in a loop.

        PLOTTING: When plotting, this graph will use the default
        spring-layout algorithm, unless a position dictionary is
        specified.

        EXAMPLES: Add one vertex to an empty graph and then show::

            sage: empty1 = graphs.EmptyGraph()
            sage: empty1.add_vertex()
            sage: empty1.show() # long time

        Use for loops to build a graph from an empty graph::

            sage: empty2 = graphs.EmptyGraph()
            sage: for i in range(5):
            ...    empty2.add_vertex() # add 5 nodes, labeled 0-4
            ...
            sage: for i in range(3):
            ...    empty2.add_edge(i,i+1) # add edges {[0:1],[1:2],[2:3]}
            ...
            sage: for i in range(4)[1:]:
            ...    empty2.add_edge(4,i) # add edges {[1:4],[2:4],[3:4]}
            ...
            sage: empty2.show() # long time
        """
        return graph.Graph(sparse=True)

    def ToroidalGrid2dGraph(self,n1,n2):
        r"""
        Returns a toroidal 2-dimensional grid graph with `n_1n_2` nodes
        (`n_1` rows and `n_2` columns).

        The toroidal 2-dimensional grid with parameters `n_1,n_2` is
        the 2-dimensional grid graph with identical parameters
        to which are added the edges `((i,0),(i,n_2-1))` and
        `((0,i),(n_1-1,i))`.

        EXAMPLE:

        The toroidal 2-dimensional grid is a regular graph, while
        the usual 2-dimensional grid is not ::

            sage: tgrid = graphs.ToroidalGrid2dGraph(8,9)
            sage: print tgrid
            Toroidal 2D Grid Graph with parameters 8,9
            sage: grid = graphs.Grid2dGraph(8,9)
            sage: grid.is_regular()
            False
            sage: tgrid.is_regular()
            True
        """

        g = self.Grid2dGraph(n1,n2)

        g.add_edges([((i,0),(i,n2-1)) for i in range(n1)] + [((0,i),(n1-1,i)) for i in range(n2)])

        g.name("Toroidal 2D Grid Graph with parameters "+str(n1)+","+str(n2))

        return g

    def Grid2dGraph(self, n1, n2):
        r"""
        Returns a `2`-dimensional grid graph with `n_1n_2` nodes (`n_1` rows and
        `n_2` columns).

        A 2d grid graph resembles a `2` dimensional grid. All inner nodes are
        connected to their `4` neighbors. Outer (non-corner) nodes are
        connected to their `3` neighbors. Corner nodes are connected to their
        2 neighbors.

        This constructor depends on NetworkX numeric labels.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, nodes are
        labelled in (row, column) pairs with `(0, 0)` in the top left corner.
        Edges will always be horizontal and vertical - another advantage of
        filling the position dictionary.

        EXAMPLES: Construct and show a grid 2d graph Rows = `5`, Columns = `7`

        ::

            sage: g = graphs.Grid2dGraph(5,7)
            sage: g.show() # long time
        """
        pos_dict = {}
        for i in range(n1):
            y = -i
            for j in range(n2):
                x = j
                pos_dict[i,j] = (x,y)
        import networkx
        G = networkx.grid_2d_graph(n1,n2)
        return graph.Graph(G, pos=pos_dict, name="2D Grid Graph")

    def GridGraph(self, dim_list):
        """
        Returns an n-dimensional grid graph.

        INPUT:


        -  ``dim_list`` - a list of integers representing the
           number of nodes to extend in each dimension.


        PLOTTING: When plotting, this graph will use the default
        spring-layout algorithm, unless a position dictionary is
        specified.

        EXAMPLES::

            sage: G = graphs.GridGraph([2,3,4])
            sage: G.show()  # long time

        ::

            sage: C = graphs.CubeGraph(4)
            sage: G = graphs.GridGraph([2,2,2,2])
            sage: C.show()  # long time
            sage: G.show()  # long time
        """
        import networkx
        dim = [int(a) for a in dim_list]
        G = networkx.grid_graph(dim)
        return graph.Graph(G, name="Grid Graph for %s"%dim)

    def HanoiTowerGraph(self, pegs, disks, labels=True, positions=True):
        r"""
        Returns the graph whose vertices are the states of the
        Tower of Hanoi puzzle, with edges representing legal moves between states.

        INPUT:

        - ``pegs`` - the number of pegs in the puzzle, 2 or greater
        - ``disks`` - the number of disks in the puzzle, 1 or greater
        - ``labels`` - default: ``True``, if ``True`` the graph contains
          more meaningful labels, see explanation below.  For large instances,
          turn off labels for much faster creation of the graph.
        - ``positions`` - default: ``True``, if ``True`` the graph contains
          layout information.  This creates a planar layout for the case
          of three pegs.  For large instances, turn off layout information
          for much faster creation of the graph.

        OUTPUT:

        The Tower of Hanoi puzzle has a certain number of identical pegs
        and a certain number of disks, each of a different radius.
        Initially the disks are all on a single peg, arranged
        in order of their radii, with the largest on the bottom.

        The goal of the puzzle is to move the disks to any other peg,
        arranged in the same order.  The one constraint is that the
        disks resident on any one peg must always be arranged with larger
        radii lower down.

        The vertices of this graph represent all the possible states
        of this puzzle.  Each state of the puzzle is a tuple with length
        equal to the number of disks, ordered by largest disk first.
        The entry of the tuple is the peg where that disk resides.
        Since disks on a given peg must go down in size as we go
        up the peg, this totally describes the state of the puzzle.

        For example ``(2,0,0)`` means the large disk is on peg 2, the
        medium disk is on peg 0, and the small disk is on peg 0
        (and we know the small disk must be above the medium disk).
        We encode these tuples as integers with a base equal to
        the number of pegs, and low-order digits to the right.

        Two vertices are adjacent if we can change the puzzle from
        one state to the other by moving a single disk.  For example,
        ``(2,0,0)`` is adjacent to ``(2,0,1)`` since we can move
        the small disk off peg 0 and onto (the empty) peg 1.
        So the solution to a 3-disk puzzle (with at least
        two pegs) can be expressed by the shortest path between
        ``(0,0,0)`` and ``(1,1,1)``.

        For greatest speed we create graphs with integer vertices,
        where we encode the tuples as integers with a base equal
        to the number of pegs, and low-order digits to the right.
        So for example, in a 3-peg puzzle with 5 disks, the
        state ``(1,2,0,1,1)`` is encoded as
        `1\ast 3^4 + 2\ast 3^3 + 0\ast 3^2 + 1\ast 3^1 + 1\ast 3^0 = 139`.

        For smaller graphs, the labels that are the tuples are informative,
        but slow down creation of the graph.  Likewise computing layout
        information also incurs a significant speed penalty. For maximum
        speed, turn off labels and layout and decode the
        vertices explicitly as needed.  The
        :meth:`sage.rings.integer.Integer.digits`
        with the ``padsto`` option is a quick way to do this, though you
        may want to reverse the list that is output.

        PLOTTING:

        The layout computed when ``positions = True`` will
        look especially good for the three-peg case, when the graph is known
        to be planar.  Except for two small cases on 4 pegs, the graph is
        otherwise not planar, and likely there is a better way to layout
        the vertices.

        EXAMPLES:

        A classic puzzle uses 3 pegs.  We solve the 5 disk puzzle using
        integer labels and report the minimum number of moves required.
        Note that `3^5-1` is the state where all 5 disks
        are on peg 2. ::

            sage: H = graphs.HanoiTowerGraph(3, 5, labels=False, positions=False)
            sage: H.distance(0, 3^5-1)
            31

        A slightly larger instance. ::

            sage: H = graphs.HanoiTowerGraph(4, 6, labels=False, positions=False)
            sage: H.num_verts()
            4096
            sage: H.distance(0, 4^6-1)
            17

        For a small graph, labels and layout information can be useful.
        Here we explicitly list a solution as a list of states. ::

            sage: H = graphs.HanoiTowerGraph(3, 3, labels=True, positions=True)
            sage: H.shortest_path((0,0,0), (1,1,1))
            [(0, 0, 0), (0, 0, 1), (0, 2, 1), (0, 2, 2), (1, 2, 2), (1, 2, 0), (1, 1, 0), (1, 1, 1)]

        Some facts about this graph with `p` pegs and `d` disks:

        - only automorphisms are the "obvious" ones - renumber the pegs.
        - chromatic number is less than or equal to `p`
        - independence number is `p^{d-1}`

        ::

            sage: H = graphs.HanoiTowerGraph(3,4,labels=False,positions=False)
            sage: H.automorphism_group().is_isomorphic(SymmetricGroup(3))
            True
            sage: H.chromatic_number()
            3
            sage: len(H.independent_set()) == 3^(4-1)
            True

        TESTS:

        It is an error to have just one peg (or less). ::

            sage: graphs.HanoiTowerGraph(1, 5)
            Traceback (most recent call last):
            ...
            ValueError: Pegs for Tower of Hanoi graph should be two or greater (not 1)

        It is an error to have zero disks (or less). ::

            sage: graphs.HanoiTowerGraph(2, 0)
            Traceback (most recent call last):
            ...
            ValueError: Disks for Tower of Hanoi graph should be one or greater (not 0)

        .. rubric:: Citations

        .. [ARETT-DOREE] Danielle Arett and Su Doree, Coloring and counting on the tower of Hanoi graphs,
           Mathematics Magazine, to appear.


        AUTHOR:

        - Rob Beezer, (2009-12-26), with assistance from Su Doree

        """

        # sanitize input
        from sage.rings.all import Integer
        pegs = Integer(pegs)
        if pegs < 2:
            raise ValueError("Pegs for Tower of Hanoi graph should be two or greater (not %d)" % pegs)
        disks = Integer(disks)
        if disks < 1:
            raise ValueError("Disks for Tower of Hanoi graph should be one or greater (not %d)" % disks)

        # Each state of the puzzle is a tuple with length
        # equal to the number of disks, ordered by largest disk first
        # The entry of the tuple is the peg where that disk resides
        # Since disks on a given peg must go down in size as we go
        # up the peg, this totally describes the puzzle
        # We encode these tuples as integers with a base equal to
        # the number of pegs, and low-order digits to the right

        # complete graph on number of pegs when just a single disk
        edges = [[i,j] for i in range(pegs) for j in range(i+1,pegs)]

        nverts = 1
        for d in range(2, disks+1):
            prevedges = edges      # remember subgraph to build from
            nverts = pegs*nverts   # pegs^(d-1)
            edges = []

            # Take an edge, change its two states in the same way by adding
            # a large disk to the bottom of the same peg in each state
            # This is accomplished by adding a multiple of pegs^(d-1)
            for p in range(pegs):
                largedisk = p*nverts
                for anedge in prevedges:
                    edges.append([anedge[0]+largedisk, anedge[1]+largedisk])

            # Two new states may only differ in the large disk
            # being the only disk on two different pegs, thus
            # otherwise being a common state with one less disk
            # We construct all such pairs of new states and add as edges
            from sage.combinat.subset import Subsets
            for state in range(nverts):
                emptypegs = range(pegs)
                reduced_state = state
                for i in range(d-1):
                    apeg = reduced_state % pegs
                    if apeg in emptypegs:
                        emptypegs.remove(apeg)
                    reduced_state = reduced_state//pegs
                for freea, freeb in Subsets(emptypegs, 2):
                    edges.append([freea*nverts+state,freeb*nverts+state])

        H = graph.Graph({}, loops=False, multiedge=False)
        H.add_edges(edges)


        # Making labels and/or computing positions can take a long time,
        # relative to just constructing the edges on integer vertices.
        # We try to minimize coercion overhead, but need Sage
        # Integers in order to use digits() for labels.
        # Getting the digits with custom code was no faster.
        # Layouts are circular (symmetric on the number of pegs)
        # radiating outward to the number of disks (radius)
        # Algorithm uses some combination of alternate
        # clockwise/counterclockwise placements, which
        # works well for three pegs (planar layout)
        #
        from sage.functions.trig import sin, cos, csc
        if labels or positions:
            mapping = {}
            pos = {}
            a = Integer(-1)
            one = Integer(1)
            if positions:
                radius_multiplier = 1 + csc(pi/pegs)
                sine = []; cosine = []
                for i in range(pegs):
                    angle = 2*i*pi/float(pegs)
                    sine.append(sin(angle))
                    cosine.append(cos(angle))
            for i in range(pegs**disks):
                a += one
                state = a.digits(base=pegs, padto=disks)
                if labels:
                    state.reverse()
                    mapping[i] = tuple(state)
                    state.reverse()
                if positions:
                    locx = 0.0; locy = 0.0
                    radius = 1.0
                    parity = -1.0
                    for index in range(disks):
                        p = state[index]
                        radius *= radius_multiplier
                        parity *= -1.0
                        locx_temp = cosine[p]*locx - parity*sine[p]*locy + radius*cosine[p]
                        locy_temp = parity*sine[p]*locx + cosine[p]*locy - radius*parity*sine[p]
                        locx = locx_temp
                        locy = locy_temp
                    pos[i] = (locx,locy)
            # set positions, then relabel (not vice versa)
            if positions:
                H.set_pos(pos)
            if labels:
                H.relabel(mapping)

        return H

    def HouseGraph(self):
        """
        Returns a house graph with 5 nodes.

        A house graph is named for its shape. It is a triangle (roof) over a
        square (walls).

        This constructor depends on NetworkX numeric labeling.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the house
        graph is drawn with the first node in the lower-left corner of the
        house, the second in the lower-right corner of the house. The third
        node is in the upper-left corner connecting the roof to the wall,
        and the fourth is in the upper-right corner connecting the roof to
        the wall. The fifth node is the top of the roof, connected only to
        the third and fourth.

        EXAMPLES: Construct and show a house graph

        ::

            sage: g = graphs.HouseGraph()
            sage: g.show() # long time
        """
        pos_dict = {0:(-1,0),1:(1,0),2:(-1,1),3:(1,1),4:(0,2)}
        import networkx
        G = networkx.house_graph()
        return graph.Graph(G, pos=pos_dict, name="House Graph")

    def HouseXGraph(self):
        """
        Returns a house X graph with 5 nodes.

        A house X graph is a house graph with two additional edges. The
        upper-right corner is connected to the lower-left. And the
        upper-left corner is connected to the lower-right.

        This constructor depends on NetworkX numeric labeling.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the house X
        graph is drawn with the first node in the lower-left corner of the
        house, the second in the lower-right corner of the house. The third
        node is in the upper-left corner connecting the roof to the wall,
        and the fourth is in the upper-right corner connecting the roof to
        the wall. The fifth node is the top of the roof, connected only to
        the third and fourth.

        EXAMPLES: Construct and show a house X graph

        ::

            sage: g = graphs.HouseXGraph()
            sage: g.show() # long time
        """
        pos_dict = {0:(-1,0),1:(1,0),2:(-1,1),3:(1,1),4:(0,2)}
        import networkx
        G = networkx.house_x_graph()
        return graph.Graph(G, pos=pos_dict, name="House Graph")

    def KrackhardtKiteGraph(self):
        """
        Returns a Krackhardt kite graph with 10 nodes.

        The Krackhardt kite graph was originally developed by David
        Krackhardt for the purpose of studying social networks. It is used
        to show the distinction between: degree centrality, betweeness
        centrality, and closeness centrality. For more information read the
        plotting section below in conjunction with the example.

        REFERENCES:

        - [1] Kreps, V. (2002). "Social Network Analysis".  [Online] Available:
          http://www.fsu.edu/~spap/water/network/intro.htm [2007,
          January 17]

        This constructor depends on NetworkX numeric labeling.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the graph is
        drawn left to right, in top to bottom row sequence of [2, 3, 2, 1,
        1, 1] nodes on each row. This places the fourth node (3) in the
        center of the kite, with the highest degree. But the fourth node
        only connects nodes that are otherwise connected, or those in its
        clique (i.e.: Degree Centrality). The eighth (7) node is where the
        kite meets the tail. It has degree = 3, less than the average, but
        is the only connection between the kite and tail (i.e.: Betweenness
        Centrality). The sixth and seventh nodes (5 and 6) are drawn in the
        third row and have degree = 5. These nodes have the shortest path
        to all other nodes in the graph (i.e.: Closeness Centrality).
        Please execute the example for visualization.

        EXAMPLE: Construct and show a Krackhardt kite graph

        ::

            sage: g = graphs.KrackhardtKiteGraph()
            sage: g.show() # long time
        """
        pos_dict = {0:(-1,4),1:(1,4),2:(-2,3),3:(0,3),4:(2,3),5:(-1,2),6:(1,2),7:(0,1),8:(0,0),9:(0,-1)}
        import networkx
        G = networkx.krackhardt_kite_graph()
        return graph.Graph(G, pos=pos_dict, name="Krackhardt Kite Graph")

    def LadderGraph(self, n):
        """
        Returns a ladder graph with 2\*n nodes.

        A ladder graph is a basic structure that is typically displayed as
        a ladder, i.e.: two parallel path graphs connected at each
        corresponding node pair.

        This constructor depends on NetworkX numeric labels.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, each ladder
        graph will be displayed horizontally, with the first n nodes
        displayed left to right on the top horizontal line.

        EXAMPLES: Construct and show a ladder graph with 14 nodes

        ::

            sage: g = graphs.LadderGraph(7)
            sage: g.show() # long time

        Create several ladder graphs in a Sage graphics array

        ::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.LadderGraph(i+2)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time
        """
        pos_dict = {}
        for i in range(n):
            pos_dict[i] = (i,1)
        for i in range(n,2*n):
            x = i - n
            pos_dict[i] = (x,0)
        import networkx
        G = networkx.ladder_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Ladder graph")

    def LollipopGraph(self, n1, n2):
        """
        Returns a lollipop graph with n1+n2 nodes.

        A lollipop graph is a path graph (order n2) connected to a complete
        graph (order n1). (A barbell graph minus one of the bells).

        This constructor depends on NetworkX numeric labels.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the complete
        graph will be drawn in the lower-left corner with the (n1)th node
        at a 45 degree angle above the right horizontal center of the
        complete graph, leading directly into the path graph.

        EXAMPLES: Construct and show a lollipop graph Candy = 13, Stick =
        4

        ::

            sage: g = graphs.LollipopGraph(13,4)
            sage: g.show() # long time

        Create several lollipop graphs in a Sage graphics array

        ::

            sage: g = []
            sage: j = []
            sage: for i in range(6):
            ...    k = graphs.LollipopGraph(i+3,4)
            ...    g.append(k)
            ...
            sage: for i in range(2):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time
        """
        pos_dict = {}

        for i in range(n1):
            x = float(cos((pi/4) - ((2*pi)/n1)*i) - n2/2 - 1)
            y = float(sin((pi/4) - ((2*pi)/n1)*i) - n2/2 - 1)
            j = n1-1-i
            pos_dict[j] = (x,y)
        for i in range(n1, n1+n2):
            x = float(i - n1 - n2/2 + 1)
            y = float(i - n1 - n2/2 + 1)
            pos_dict[i] = (x,y)

        import networkx
        G = networkx.lollipop_graph(n1,n2)
        return graph.Graph(G, pos=pos_dict, name="Lollipop Graph")

    def PathGraph(self, n, pos=None):
        """
        Returns a path graph with n nodes. Pos argument takes a string
        which is either 'circle' or 'line', (otherwise the default is
        used). See the plotting section below for more detail.

        A path graph is a graph where all inner nodes are connected to
        their two neighbors and the two end-nodes are connected to their
        one inner neighbors. (i.e.: a cycle graph without the first and
        last node connected).

        This constructor depends on NetworkX numeric labels.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the graph may
        be drawn in one of two ways: The 'line' argument will draw the
        graph in a horizontal line (left to right) if there are less than
        11 nodes. Otherwise the 'line' argument will append horizontal
        lines of length 10 nodes below, alternating left to right and right
        to left. The 'circle' argument will cause the graph to be drawn in
        a cycle-shape, with the first node at the top and then about the
        circle in a clockwise manner. By default (without an appropriate
        string argument) the graph will be drawn as a 'circle' if 10 n 41
        and as a 'line' for all other n.

        EXAMPLES: Show default drawing by size: 'line': n 11

        ::

            sage: p = graphs.PathGraph(10)
            sage: p.show() # long time

        'circle': 10 n 41

        ::

            sage: q = graphs.PathGraph(25)
            sage: q.show() # long time

        'line': n 40

        ::

            sage: r = graphs.PathGraph(55)
            sage: r.show() # long time

        Override the default drawing::

            sage: s = graphs.PathGraph(5,'circle')
            sage: s.show() # long time
        """
        pos_dict = {}

        # Choose appropriate drawing pattern
        circle = False
        if pos == "circle": circle = True
        elif pos == "line": circle = False
        # Otherwise use default by size of n
        elif 10 < n < 41: circle = True

        # Draw 'circle'
        if circle:
            for i in range(n):
                x = float(cos((pi/2) + ((2*pi)/n)*i))
                y = float(sin((pi/2) + ((2*pi)/n)*i))
                pos_dict[i] = (x,y)
        # Draw 'line'
        else:
            counter = 0 # node index
            rem = n%10 # remainder to appear on last row
            rows = n//10 # number of rows (not counting last row)
            lr = True # left to right

            for i in range(rows): # note that rows doesn't include last row
                y = -i
                for j in range(10):
                    if lr:
                        x = j
                    else:
                        x = 9 - j
                    pos_dict[counter] = (x,y)
                    counter += 1
                if lr: lr = False
                else: lr = True
            y = -rows
            for j in range(rem): # last row
                if lr:
                    x = j
                else:
                    x = 9 - j
                pos_dict[counter] = (x,y)
                counter += 1

        import networkx
        G = networkx.path_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Path Graph")

    def StarGraph(self, n):
        """
        Returns a star graph with n+1 nodes.

        A Star graph is a basic structure where one node is connected to
        all other nodes.

        This constructor is dependent on NetworkX numeric labels.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, each star
        graph will be displayed with the first (0) node in the center, the
        second node (1) at the top, with the rest following in a
        counterclockwise manner. (0) is the node connected to all other
        nodes.

        The star graph is a good opportunity to compare efficiency of
        filling a position dictionary vs. using the spring-layout algorithm
        for plotting. As far as display, the spring-layout should push all
        other nodes away from the (0) node, and thus look very similar to
        this constructor's positioning.

        EXAMPLES::

            sage: import networkx

        Compare the plots::

            sage: n = networkx.star_graph(23)
            sage: spring23 = Graph(n)
            sage: posdict23 = graphs.StarGraph(23)
            sage: spring23.show() # long time
            sage: posdict23.show() # long time

        View many star graphs as a Sage Graphics Array

        With this constructor (i.e., the position dictionary filled)

        ::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.StarGraph(i+3)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time

        Compared to plotting with the spring-layout algorithm

        ::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    spr = networkx.star_graph(i+3)
            ...    k = Graph(spr)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time
        """
        pos_dict = {}
        pos_dict[0] = (0,0)
        for i in range(1,n+1):
            x = float(cos((pi/2) + ((2*pi)/n)*(i-1)))
            y = float(sin((pi/2) + ((2*pi)/n)*(i-1)))
            pos_dict[i] = (x,y)
        import networkx
        G = networkx.star_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Star graph")

    def WheelGraph(self, n):
        """
        Returns a Wheel graph with n nodes.

        A Wheel graph is a basic structure where one node is connected to
        all other nodes and those (outer) nodes are connected cyclically.

        This constructor depends on NetworkX numeric labels.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, each wheel
        graph will be displayed with the first (0) node in the center, the
        second node at the top, and the rest following in a
        counterclockwise manner.

        With the wheel graph, we see that it doesn't take a very large n at
        all for the spring-layout to give a counter-intuitive display. (See
        Graphics Array examples below).

        EXAMPLES: We view many wheel graphs with a Sage Graphics Array,
        first with this constructor (i.e., the position dictionary
        filled)::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.WheelGraph(i+3)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time

        Next, using the spring-layout algorithm::

            sage: import networkx
            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    spr = networkx.wheel_graph(i+3)
            ...    k = Graph(spr)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time

        Compare the plotting::

            sage: n = networkx.wheel_graph(23)
            sage: spring23 = Graph(n)
            sage: posdict23 = graphs.WheelGraph(23)
            sage: spring23.show() # long time
            sage: posdict23.show() # long time
        """
        pos_dict = {}
        pos_dict[0] = (0,0)
        for i in range(1,n):
            x = float(cos((pi/2) + ((2*pi)/(n-1))*(i-1)))
            y = float(sin((pi/2) + ((2*pi)/(n-1))*(i-1)))
            pos_dict[i] = (x,y)
        import networkx
        G = networkx.wheel_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Wheel graph")

################################################################################
#   Platonic Solids
################################################################################

    def TetrahedralGraph(self):
        """
        Returns a tetrahedral graph (with 4 nodes).

        A tetrahedron is a 4-sided triangular pyramid. The tetrahedral
        graph corresponds to the connectivity of the vertices of the
        tetrahedron. This graph is equivalent to a wheel graph with 4 nodes
        and also a complete graph on four nodes. (See examples below).

        PLOTTING: The tetrahedral graph should be viewed in 3 dimensions.
        We chose to use the default spring-layout algorithm here, so that
        multiple iterations might yield a different point of reference for
        the user. We hope to add rotatable, 3-dimensional viewing in the
        future. In such a case, a string argument will be added to select
        the flat spring-layout over a future implementation.

        EXAMPLES: Construct and show a Tetrahedral graph

        ::

            sage: g = graphs.TetrahedralGraph()
            sage: g.show() # long time

        The following example requires networkx::

            sage: import networkx as NX

        Compare this Tetrahedral, Wheel(4), Complete(4), and the
        Tetrahedral plotted with the spring-layout algorithm below in a
        Sage graphics array::

            sage: tetra_pos = graphs.TetrahedralGraph()
            sage: tetra_spring = Graph(NX.tetrahedral_graph())
            sage: wheel = graphs.WheelGraph(4)
            sage: complete = graphs.CompleteGraph(4)
            sage: g = [tetra_pos, tetra_spring, wheel, complete]
            sage: j = []
            sage: for i in range(2):
            ...    n = []
            ...    for m in range(2):
            ...        n.append(g[i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time
        """
        import networkx
        G = networkx.tetrahedral_graph()
        return graph.Graph(G, name="Tetrahedron")

    def HexahedralGraph(self):
        """
        Returns a hexahedral graph (with 8 nodes).

        A regular hexahedron is a 6-sided cube. The hexahedral graph
        corresponds to the connectivity of the vertices of the hexahedron.
        This graph is equivalent to a 3-cube.

        PLOTTING: The hexahedral graph should be viewed in 3 dimensions. We
        chose to use the default spring-layout algorithm here, so that
        multiple iterations might yield a different point of reference for
        the user. We hope to add rotatable, 3-dimensional viewing in the
        future. In such a case, a string argument will be added to select
        the flat spring-layout over a future implementation.

        EXAMPLES: Construct and show a Hexahedral graph

        ::

            sage: g = graphs.HexahedralGraph()
            sage: g.show() # long time

        Create several hexahedral graphs in a Sage graphics array. They
        will be drawn differently due to the use of the spring-layout
        algorithm.

        ::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.HexahedralGraph()
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time
        """
        return graph.Graph({0:[1,3,4], 1:[2,5], 2:[3,6], 3:[7], 4:[5,7],\
                            5:[6], 6:[7]}, name="Hexahedron")

    def OctahedralGraph(self):
        """
        Returns an Octahedral graph (with 6 nodes).

        The regular octahedron is an 8-sided polyhedron with triangular
        faces. The octahedral graph corresponds to the connectivity of the
        vertices of the octahedron. It is the line graph of the tetrahedral
        graph. The octahedral is symmetric, so the spring-layout algorithm
        will be very effective for display.

        PLOTTING: The Octahedral graph should be viewed in 3 dimensions. We
        chose to use the default spring-layout algorithm here, so that
        multiple iterations might yield a different point of reference for
        the user. We hope to add rotatable, 3-dimensional viewing in the
        future. In such a case, a string argument will be added to select
        the flat spring-layout over a future implementation.

        EXAMPLES: Construct and show an Octahedral graph

        ::

            sage: g = graphs.OctahedralGraph()
            sage: g.show() # long time

        Create several octahedral graphs in a Sage graphics array They will
        be drawn differently due to the use of the spring-layout algorithm

        ::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.OctahedralGraph()
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time
        """
        import networkx
        G = networkx.octahedral_graph()
        return graph.Graph(G, name="Octahedron")

    def IcosahedralGraph(self):
        """
        Returns an Icosahedral graph (with 12 nodes).

        The regular icosahedron is a 20-sided triangular polyhedron. The
        icosahedral graph corresponds to the connectivity of the vertices
        of the icosahedron. It is dual to the dodecahedral graph. The
        icosahedron is symmetric, so the spring-layout algorithm will be
        very effective for display.

        PLOTTING: The Icosahedral graph should be viewed in 3 dimensions.
        We chose to use the default spring-layout algorithm here, so that
        multiple iterations might yield a different point of reference for
        the user. We hope to add rotatable, 3-dimensional viewing in the
        future. In such a case, a string argument will be added to select
        the flat spring-layout over a future implementation.

        EXAMPLES: Construct and show an Octahedral graph

        ::

            sage: g = graphs.IcosahedralGraph()
            sage: g.show() # long time

        Create several icosahedral graphs in a Sage graphics array. They
        will be drawn differently due to the use of the spring-layout
        algorithm.

        ::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.IcosahedralGraph()
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time
        """
        import networkx
        G = networkx.icosahedral_graph()
        return graph.Graph(G, name="Icosahedron")

    def DodecahedralGraph(self):
        """
        Returns a Dodecahedral graph (with 20 nodes)

        The dodecahedral graph is cubic symmetric, so the spring-layout
        algorithm will be very effective for display. It is dual to the
        icosahedral graph.

        PLOTTING: The Dodecahedral graph should be viewed in 3 dimensions.
        We chose to use the default spring-layout algorithm here, so that
        multiple iterations might yield a different point of reference for
        the user. We hope to add rotatable, 3-dimensional viewing in the
        future. In such a case, a string argument will be added to select
        the flat spring-layout over a future implementation.

        EXAMPLES: Construct and show a Dodecahedral graph

        ::

            sage: g = graphs.DodecahedralGraph()
            sage: g.show() # long time

        Create several dodecahedral graphs in a Sage graphics array They
        will be drawn differently due to the use of the spring-layout
        algorithm

        ::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.DodecahedralGraph()
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time
        """
        import networkx
        G = networkx.dodecahedral_graph()
        return graph.Graph(G, name="Dodecahedron")

    #######################################################################
    #   Named Graphs
    #######################################################################

    def ChvatalGraph(self):
        r"""
        Returns the Chvatal graph.

        Chvatal graph is one of the few known graphs to satisfy Grunbaum's
        conjecture that for every m, n, there is an m-regular,
        m-chromatic graph of girth at least n. For more information, see this
        `Wikipedia article on the Chvatal graph <http://en.wikipedia.org/wiki/Chv%C3%A1tal_graph>`_.

        EXAMPLES:

        The Chvatal graph has 12 vertices and 24 edges. It is a 4-regular,
        4-chromatic graph with radius 2, diameter 2, and girth 4. ::

            sage: G = graphs.ChvatalGraph(); G
            Chvatal graph: Graph on 12 vertices
            sage: G.order(); G.size()
            12
            24
            sage: G.degree()
            [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
            sage: G.chromatic_number()
            4
            sage: G.radius(); G.diameter(); G.girth()
            2
            2
            4
        """
        import networkx
        pos_dict = {}
        for i in range(5, 10):
            x = float(cos((pi / 2) + ((2 * pi) / 5) * i))
            y = float(sin((pi / 2) + ((2 * pi) / 5) * i))
            pos_dict[i] = (x, y)
        for i in range(5):
            x = float(2 * (cos((pi / 2) + ((2 * pi) / 5) * (i - 5))))
            y = float(2 * (sin((pi / 2) + ((2 * pi) / 5) * (i - 5))))
            pos_dict[i] = (x, y)
        pos_dict[10] = (0.5, 0)
        pos_dict[11] = (-0.5, 0)

        return graph.Graph(networkx.chvatal_graph(), pos=pos_dict, name="Chvatal graph")

    def DesarguesGraph(self):
        """
        Returns the Desargues graph.

        PLOTTING: The layout chosen is the same as on the cover of [1].

        EXAMPLE::

            sage: D = graphs.DesarguesGraph()
            sage: L = graphs.LCFGraph(20,[5,-5,9,-9],5)
            sage: D.is_isomorphic(L)
            True
            sage: D.show()  # long time

        REFERENCE:

        - [1] Harary, F. Graph Theory. Reading, MA: Addison-Wesley,
          1994.
        """
        G=graphs.GeneralizedPetersenGraph(10,3)
        G.name("Desargues Graph")
        return G

    def FlowerSnark(self):
        """
        Returns a Flower Snark.

        A flower snark has 20 vertices. It is part of the class of
        biconnected cubic graphs with edge chromatic number = 4, known as
        snarks. (i.e.: the Petersen graph). All snarks are not Hamiltonian,
        non-planar and have Petersen graph graph minors.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the nodes are
        drawn 0-14 on the outer circle, and 15-19 in an inner pentagon.

        REFERENCES:

        - [1] Weisstein, E. (1999). "Flower Snark - from Wolfram
          MathWorld". [Online] Available:
          http://mathworld.wolfram.com/FlowerSnark.html [2007, February 17]

        EXAMPLES: Inspect a flower snark::

            sage: F = graphs.FlowerSnark()
            sage: F
            Flower Snark: Graph on 20 vertices
            sage: F.graph6_string()
            'ShCGHC@?GGg@?@?Gp?K??C?CA?G?_G?Cc'

        Now show it::

            sage: F.show() # long time
        """
        pos_dict = {}
        for i in range(15):
            x = float(2.5*(cos((pi/2) + ((2*pi)/15)*i)))
            y = float(2.5*(sin((pi/2) + ((2*pi)/15)*i)))
            pos_dict[i] = (x,y)
        for i in range(15,20):
            x = float(cos((pi/2) + ((2*pi)/5)*i))
            y = float(sin((pi/2) + ((2*pi)/5)*i))
            pos_dict[i] = (x,y)
        return graph.Graph({0:[1,14,15],1:[2,11],2:[3,7],3:[2,4,16],4:[5,14], \
                            5:[6,10],6:[5,7,17],8:[7,9,13],9:[10,18],11:[10,12], \
                            12:[13,19],13:[14],15:[19],16:[15,17],18:[17,19]}, \
                            pos=pos_dict, name="Flower Snark")

    def FruchtGraph(self):
        """
        Returns a Frucht Graph.

        A Frucht graph has 12 nodes and 18 edges. It is the smallest cubic
        identity graph. It is planar and it is Hamiltonian.

        This constructor is dependent on NetworkX's numeric labeling.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the first
        seven nodes are on the outer circle, with the next four on an inner
        circle and the last in the center.

        REFERENCES:

        - [1] Weisstein, E. (1999). "Frucht Graph - from Wolfram
          MathWorld". [Online] Available:
          http://mathworld.wolfram.com/FruchtGraph.html [2007, February 17]

        EXAMPLES::

            sage: FRUCHT = graphs.FruchtGraph()
            sage: FRUCHT
            Frucht graph: Graph on 12 vertices
            sage: FRUCHT.graph6_string()
            'KhCKM?_EGK?L'
            sage: (graphs.FruchtGraph()).show() # long time
        """
        pos_dict = {}
        for i in range(7):
            x = float(2*(cos((pi/2) + ((2*pi)/7)*i)))
            y = float(2*(sin((pi/2) + ((2*pi)/7)*i)))
            pos_dict[i] = (x,y)
        pos_dict[7] = (0,1)
        pos_dict[8] = (-1,0)
        pos_dict[9] = (0,-1)
        pos_dict[10] = (1,0)
        pos_dict[11] = (0,0)
        import networkx
        G = networkx.frucht_graph()
        return graph.Graph(G, pos=pos_dict, name="Frucht graph")

    def HeawoodGraph(self):
        """
        Returns a Heawood graph.

        The Heawood graph is a cage graph that has 14 nodes. It is a cubic
        symmetric graph. (See also the Moebius-Kantor graph). It is
        nonplanar and Hamiltonian. It has diameter = 3, radius = 3, girth =
        6, chromatic number = 2. It is 4-transitive but not 5-transitive.

        This constructor is dependent on NetworkX's numeric labeling.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the nodes are
        positioned in a circular layout with the first node appearing at
        the top, and then continuing counterclockwise.

        REFERENCES:

        - [1] Weisstein, E. (1999). "Heawood Graph - from Wolfram
          MathWorld". [Online] Available:
          http://mathworld.wolfram.com/HeawoodGraph.html [2007, February 17]

        EXAMPLES::

            sage: H = graphs.HeawoodGraph()
            sage: H
            Heawood graph: Graph on 14 vertices
            sage: H.graph6_string()
            'MhEGHC@AI?_PC@_G_'
            sage: (graphs.HeawoodGraph()).show() # long time
        """
        pos_dict = {}
        for i in range(14):
            x = float(cos((pi/2) + (pi/7)*i))
            y = float(sin((pi/2) + (pi/7)*i))
            pos_dict[i] = (x,y)
        import networkx
        G = networkx.heawood_graph()
        return graph.Graph(G, pos=pos_dict, name="Heawood graph")

    def HigmanSimsGraph(self, relabel=True):
        r"""
        The Higman-Sims graph is a remarkable strongly regular
        graph of degree 22 on 100 vertices.  For example, it can
        be split into two sets of 50 vertices each, so that each
        half induces a subgraph isomorphic to the
        Hoffman-Singleton graph
        (:meth:`~HoffmanSingletonGraph`).
        This can be done in 352 ways (see [BROUWER-HS-2009]_).

        Its most famous property is that the automorphism
        group has an index 2 subgroup which is one of the
        26 sporadic groups. [HIGMAN1968]_

        The construction used here follows [HAFNER2004]_.

        INPUT:

        - ``relabel`` - default: ``True``.  If ``True`` the
          vertices will be labeled with consecutive integers.
          If ``False`` the labels are strings that are three
          digits long. "xyz" means the vertex is in group
          x (zero through three), pentagon or pentagram y
          (zero through four), and is vertex z (zero
          through four) of that pentagon or pentagram.
          See [HAFNER2004]_ for more.

        OUTPUT:

        The Higman-Sims graph.

        EXAMPLES:

        A split into the first 50 and last 50 vertices
        will induce two copies of the Hoffman-Singleton graph,
        and we illustrate another such split, which is obvious
        based on the construction used. ::

            sage: H = graphs.HigmanSimsGraph()
            sage: A = H.subgraph(range(0,50))
            sage: B = H.subgraph(range(50,100))
            sage: K = graphs.HoffmanSingletonGraph()
            sage: K.is_isomorphic(A) and K.is_isomorphic(B)
            True
            sage: C = H.subgraph(range(25,75))
            sage: D = H.subgraph(range(0,25)+range(75,100))
            sage: K.is_isomorphic(C) and K.is_isomorphic(D)
            True

        The automorphism group contains only one nontrivial
        proper normal subgroup, which is of index 2 and is
        simple.  It is known as the Higman-Sims group.  ::

            sage: H = graphs.HigmanSimsGraph()
            sage: G = H.automorphism_group()
            sage: g=G.order(); g
            88704000
            sage: K = G.normal_subgroups()[1]
            sage: K.is_simple()
            True
            sage: g//K.order()
            2

        REFERENCES:

            .. [BROUWER-HS-2009] `Higman-Sims graph
               <http://www.win.tue.nl/~aeb/graphs/Higman-Sims.html>`_.
               Andries E. Brouwer, accessed 24 October 2009.
            .. [HIGMAN1968] A simple group of order 44,352,000,
               Math.Z. 105 (1968) 110-113. D.G. Higman & C. Sims.
            .. [HAFNER2004] `On the graphs of Hoffman-Singleton and
               Higman-Sims
               <http://www.combinatorics.org/Volume_11/PDF/v11i1r77.pdf>`_.
               The Electronic Journal of Combinatorics 11 (2004), #R77,
               Paul R. Hafner, accessed 24 October 2009.

        AUTHOR:

            - Rob Beezer (2009-10-24)
        """
        HS = graph.Graph()
        HS.name('Higman-Sims graph')

        # Four groups of either five pentagons, or five pentagrams
        # 4 x 5 x 5 = 100 vertices
        # First digit is "group", second is "penta{gon|gram}", third is "vertex"
        vlist = ['%d%d%d'%(g,p,v)
                        for g in range(4) for p in range(5) for v in range(5)]
        for avertex in vlist:
            HS.add_vertex(avertex)

        # Edges: Within groups 0 and 2, joined as pentagons
        # Edges: Within groups 1 and 3, joined as pentagrams
        for g in range(4):
            shift = 1
            if g in [1,3]:
                shift += 1
            for p in range(5):
                for v in range(5):
                    HS.add_edge(('%d%d%d'%(g,p,v), '%d%d%d'%(g,p,(v+shift)%5)))

        # Edges: group 0 to group 1
        for x in range(5):
            for m in range(5):
                for c in range(5):
                    y = (m*x+c)%5
                    HS.add_edge(('0%d%d'%(x,y), '1%d%d'%(m,c)))

        # Edges: group 1 to group 2
        for m in range(5):
            for A in range(5):
                for B in range(5):
                    c = (2*(m-A)*(m-A)+B)%5
                    HS.add_edge(('1%d%d'%(m,c), '2%d%d'%(A,B)))

        # Edges: group 2 to group 3
        for A in range(5):
            for a in range(5):
                for b in range(5):
                    B = (2*A*A+3*a*A-a*a+b)%5
                    HS.add_edge(('2%d%d'%(A,B), '3%d%d'%(a,b)))

        # Edges: group 3 to group 0
        for a in range(5):
            for b in range(5):
                for x in range(5):
                    y = ((x-a)*(x-a)+b)%5
                    HS.add_edge(('3%d%d'%(a,b), '0%d%d'%(x,y)))

        # Edges: group 0 to group 2
        for x in range(5):
            for A in range(5):
                for B in range(5):
                    y = (3*x*x+A*x+B+1)%5
                    HS.add_edge(('0%d%d'%(x,y), '2%d%d'%(A,B)))
                    y = (3*x*x+A*x+B-1)%5
                    HS.add_edge(('0%d%d'%(x,y), '2%d%d'%(A,B)))

        # Edges: group 1 to group 3
        for m in range(5):
            for a in range(5):
                for b in range(5):
                    c = (m*(m-a)+b+2)%5
                    HS.add_edge(('1%d%d'%(m,c), '3%d%d'%(a,b)))
                    c = (m*(m-a)+b-2)%5
                    HS.add_edge(('1%d%d'%(m,c), '3%d%d'%(a,b)))

        # Rename to integer vertex labels, creating dictionary
        # Or not, and create identity mapping
        if relabel:
            vmap = HS.relabel(return_map=True)
        else:
            vmap={}
            for v in vlist:
                vmap[v] = v
        # Layout vertices in a circle
        # In the order given in vlist
        # Using labels from vmap
        pos_dict = {}
        for i in range(100):
            x = float(cos((pi/2) + ((2*pi)/100)*i))
            y = float(sin((pi/2) + ((2*pi)/100)*i))
            pos_dict[vmap[vlist[i]]] = (x,y)
        HS.set_pos(pos_dict)
        return HS

    def HoffmanSingletonGraph(self):
        r"""
        Returns the Hoffman-Singleton graph.

        The Hoffman-Singleton graph is the Moore graph of degree 7,
        diameter 2 and girth 5. The Hoffman-Singleton theorem states that
        any Moore graph with girth 5 must have degree 2, 3, 7 or 57. The
        first three respectively are the pentagon, the Petersen graph, and
        the Hoffman-Singleton graph. The existence of a Moore graph with
        girth 5 and degree 57 is still open.

        A Moore graph is a graph with diameter `d` and girth
        `2d + 1`. This implies that the graph is regular, and
        distance regular.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. A novel algorithm written by
        Tom Boothby gives a random layout which is pleasing to the eye.

        REFERENCES:

        - [1] Godsil, C. and Royle, G. Algebraic Graph Theory.
          Springer, 2001.

        EXAMPLES::

            sage: HS = graphs.HoffmanSingletonGraph()
            sage: Set(HS.degree())
            {7}
            sage: HS.girth()
            5
            sage: HS.diameter()
            2
            sage: HS.num_verts()
            50

        Note that you get a different layout each time you create the graph.

            sage: HS.layout()[1]
            (-0.844..., 0.535...)
            sage: graphs.HoffmanSingletonGraph().layout()[1]
            (-0.904..., 0.425...)
        """
        H = graph.Graph({ \
        'q00':['q01'], 'q01':['q02'], 'q02':['q03'], 'q03':['q04'], 'q04':['q00'], \
        'q10':['q11'], 'q11':['q12'], 'q12':['q13'], 'q13':['q14'], 'q14':['q10'], \
        'q20':['q21'], 'q21':['q22'], 'q22':['q23'], 'q23':['q24'], 'q24':['q20'], \
        'q30':['q31'], 'q31':['q32'], 'q32':['q33'], 'q33':['q34'], 'q34':['q30'], \
        'q40':['q41'], 'q41':['q42'], 'q42':['q43'], 'q43':['q44'], 'q44':['q40'], \
        'p00':['p02'], 'p02':['p04'], 'p04':['p01'], 'p01':['p03'], 'p03':['p00'], \
        'p10':['p12'], 'p12':['p14'], 'p14':['p11'], 'p11':['p13'], 'p13':['p10'], \
        'p20':['p22'], 'p22':['p24'], 'p24':['p21'], 'p21':['p23'], 'p23':['p20'], \
        'p30':['p32'], 'p32':['p34'], 'p34':['p31'], 'p31':['p33'], 'p33':['p30'], \
        'p40':['p42'], 'p42':['p44'], 'p44':['p41'], 'p41':['p43'], 'p43':['p40']})
        for j in range(5):
            for i in range(5):
                for k in range(5):
                    con = (i+j*k)%5
                    H.add_edge(('q%d%d'%(k,con),'p%d%d'%(j,i)))
        H.name('Hoffman-Singleton graph')
        from sage.combinat.combinat import permutations
        from sage.misc.prandom import randint
        P = permutations([1,2,3,4])
        qpp = [0]+P[randint(0,23)]
        ppp = [0]+P[randint(0,23)]
        qcycle = lambda i,s : ['q%s%s'%(i,(j+s)%5) for j in qpp]
        pcycle = lambda i,s : ['p%s%s'%(i,(j+s)%5) for j in ppp]
        l = 0
        s = 0
        D = []
        while l < 5:
            for q in qcycle(l,s):
                D.append(q)
            vv = 'p%s'%q[1]
            s = int([v[-1] for v in H.neighbors(q) if v[:2] == vv][0])
            for p in pcycle(l,s):
                D.append(p)
            vv = 'q%s'%(int(p[1])+1)
            v = [v[-1] for v in H.neighbors(p) if v[:2] == vv]
            if len(v):
                s = int(v[0])
            l+=1
        map = H.relabel(return_map=True)
        pos_dict = {}
        for i in range(50):
            x = float(cos((pi/2) + ((2*pi)/50)*i))
            y = float(sin((pi/2) + ((2*pi)/50)*i))
            pos_dict[map[D[i]]] = (x,y)
        H.set_pos(pos_dict)
        return H

    def KneserGraph(self,n,k):
        r"""
        Returns the Kneser Graph with parameters `n, k`.

        The Kneser Graph with parameters `n,k` is the graph
        whose vertices are the `k`-subsets of `[0,1,\dots,n-1]`, and such
        that two vertices are adjacent if their corresponding sets
        are disjoint.

        For example, the Petersen Graph can be defined
        as the Kneser Graph with parameters `5,2`.

        EXAMPLE::

            sage: KG=graphs.KneserGraph(5,2)
            sage: print KG.vertices()
            [{4, 5}, {1, 3}, {2, 5}, {2, 3}, {3, 4}, {3, 5}, {1, 4}, {1, 5}, {1, 2}, {2, 4}]
            sage: P=graphs.PetersenGraph()
            sage: P.is_isomorphic(KG)
            True

        TESTS::

            sage: KG=graphs.KneserGraph(0,0)
            Traceback (most recent call last):
            ...
            ValueError: Parameter n should be a strictly positive integer
            sage: KG=graphs.KneserGraph(5,6)
            Traceback (most recent call last):
            ...
            ValueError: Parameter k should be a strictly positive integer inferior to n
        """

        if not n>0:
            raise ValueError, "Parameter n should be a strictly positive integer"
        if not (k>0 and k<=n):
            raise ValueError, "Parameter k should be a strictly positive integer inferior to n"

        g=graph.Graph(name="Kneser graph with parameters "+str(n)+","+str(k))
        from sage.combinat.subset import Subsets

        if k>n/2:
            g.add_vertices(Subsets(n,k).list())

        S = Subsets(n,k)
        for s in S:
            for t in Subsets(S.s.difference(s),k):
                g.add_edge(s,t)

        return g

    def OddGraph(self,n):
        r"""
        Returns the Odd Graph with parameter `n`.

        The Odd Graph with parameter `n` is defined as the
        Kneser Graph with parameters `2n-1,n-1`.
        Equivalently, the Odd Graph is the graph whose vertices
        are the `n-1`-subsets of `[0,1,\dots,2(n-1)]`, and such
        that two vertices are adjacent if their corresponding sets
        are disjoint.

        For example, the Petersen Graph can be defined
        as the Odd Graph with parameter `3`.

        EXAMPLE::

            sage: OG=graphs.OddGraph(3)
            sage: print OG.vertices()
            [{4, 5}, {1, 3}, {2, 5}, {2, 3}, {3, 4}, {3, 5}, {1, 4}, {1, 5}, {1, 2}, {2, 4}]
            sage: P=graphs.PetersenGraph()
            sage: P.is_isomorphic(OG)
            True

        TESTS::

            sage: KG=graphs.OddGraph(1)
            Traceback (most recent call last):
            ...
            ValueError: Parameter n should be an integer strictly greater than 1
        """

        if not n>1:
            raise ValueError, "Parameter n should be an integer strictly greater than 1"
        g = self.KneserGraph(2*n-1,n-1)
        g.name("Odd Graph with parameter %s" % n)
        return g


    def MoebiusKantorGraph(self):
        """
        Returns a Moebius-Kantor Graph.

        A Moebius-Kantor graph is a cubic symmetric graph. (See also the
        Heawood graph). It has 16 nodes and 24 edges. It is nonplanar and
        Hamiltonian. It has diameter = 4, girth = 6, and chromatic number =
        2. It is identical to the Generalized Petersen graph, P[8,3].

        PLOTTING: See the plotting section for the generalized Petersen graphs.

        REFERENCES:

        - [1] Weisstein, E. (1999). "Moebius-Kantor Graph - from
          Wolfram MathWorld". [Online] Available:
          http://mathworld.wolfram.com/Moebius-KantorGraph.html [2007,
          February 17]

        EXAMPLES::

            sage: MK = graphs.MoebiusKantorGraph()
            sage: MK
            Moebius-Kantor Graph: Graph on 16 vertices
            sage: MK.graph6_string()
            'OhCGKE?O@?ACAC@I?Q_AS'
            sage: (graphs.MoebiusKantorGraph()).show() # long time
        """
        G=graphs.GeneralizedPetersenGraph(8,3)
        G.name("Moebius-Kantor Graph")
        return G

    def PappusGraph(self):
        """
        Returns the Pappus graph, a graph on 18 vertices.

        The Pappus graph is cubic, symmetric, and distance-regular.

        EXAMPLES::

            sage: G = graphs.PappusGraph()
            sage: G.show()  # long time
            sage: L = graphs.LCFGraph(18, [5,7,-7,7,-7,-5], 3)
            sage: L.show()  # long time
            sage: G.is_isomorphic(L)
            True
        """
        pos_dict = {}
        for i in range(6):
            pos_dict[i] = [float(cos(pi/2 + ((2*pi)/6)*i)),\
                           float(sin(pi/2 + ((2*pi)/6)*i))]
            pos_dict[6 + i] = [(2/3.0)*float(cos(pi/2 + ((2*pi)/6)*i)),\
                               (2/3.0)*float(sin(pi/2 + ((2*pi)/6)*i))]
            pos_dict[12 + i] = [(1/3.0)*float(cos(pi/2 + ((2*pi)/6)*i)),\
                                (1/3.0)*float(sin(pi/2 + ((2*pi)/6)*i))]
        return graph.Graph({0:[1,5,6],1:[2,7],2:[3,8],3:[4,9],4:[5,10],\
                            5:[11],6:[13,17],7:[12,14],8:[13,15],9:[14,16],\
                            10:[15,17],11:[12,16],12:[15],13:[16],14:[17]},\
                           pos=pos_dict, name="Pappus Graph")

    def PetersenGraph(self):
        """
        The Petersen Graph is a named graph that consists of 10 vertices
        and 15 edges, usually drawn as a five-point star embedded in a
        pentagon.

        The Petersen Graph is a common counterexample. For example, it is
        not Hamiltonian.

        PLOTTING: See the plotting section for the generalized Petersen graphs.

        EXAMPLES: We compare below the Petersen graph with the default
        spring-layout versus a planned position dictionary of [x,y]
        tuples::

            sage: petersen_spring = Graph({0:[1,4,5], 1:[0,2,6], 2:[1,3,7], 3:[2,4,8], 4:[0,3,9], 5:[0,7,8], 6:[1,8,9], 7:[2,5,9], 8:[3,5,6], 9:[4,6,7]})
            sage: petersen_spring.show() # long time
            sage: petersen_database = graphs.PetersenGraph()
            sage: petersen_database.show() # long time
        """
        P=graphs.GeneralizedPetersenGraph(5,2)
        P.name("Petersen graph")
        return P

    def ThomsenGraph(self):
        """
        Returns the Thomsen Graph.

        The Thomsen Graph is actually a complete bipartite graph with (n1,
        n2) = (3, 3). It is also called the Utility graph.

        PLOTTING: See CompleteBipartiteGraph.

        EXAMPLES::

            sage: T = graphs.ThomsenGraph()
            sage: T
            Thomsen graph: Graph on 6 vertices
            sage: T.graph6_string()
            'EFz_'
            sage: (graphs.ThomsenGraph()).show() # long time
        """
        pos_dict = {0:(-1,1),1:(0,1),2:(1,1),3:(-1,0),4:(0,0),5:(1,0)}
        import networkx
        G = networkx.complete_bipartite_graph(3,3)
        return graph.Graph(G, pos=pos_dict, name="Thomsen graph")

################################################################################
#   Families of Graphs
################################################################################

    def CirculantGraph(self, n, adjacency):
        r"""
        Returns a circulant graph with n nodes.

        A circulant graph has the property that the vertex i is connected
        with the vertices i+j and i-j for each j in adj.

        INPUT:


        -  ``n`` - number of vertices in the graph

        -  ``adjacency`` - the list of j values


        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, each circulant
        graph will be displayed with the first (0) node at the top, with
        the rest following in a counterclockwise manner.

        Filling the position dictionary in advance adds O(n) to the
        constructor.

        EXAMPLES: Compare plotting using the predefined layout and
        networkx::

            sage: import networkx
            sage: n = networkx.cycle_graph(23)
            sage: spring23 = Graph(n)
            sage: posdict23 = graphs.CirculantGraph(23,2)
            sage: spring23.show() # long time
            sage: posdict23.show() # long time

        We next view many cycle graphs as a Sage graphics array. First we
        use the ``CirculantGraph`` constructor, which fills in
        the position dictionary::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.CirculantGraph(i+3,i)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time

        Compare to plotting with the spring-layout algorithm::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    spr = networkx.cycle_graph(i+3)
            ...    k = Graph(spr)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time

        Passing a 1 into adjacency should give the cycle.

        ::

            sage: graphs.CirculantGraph(6,1)==graphs.CycleGraph(6)
            True
            sage: graphs.CirculantGraph(7,[1,3]).edges(labels=false)
            [(0, 1),
            (0, 3),
            (0, 4),
            (0, 6),
            (1, 2),
            (1, 4),
            (1, 5),
            (2, 3),
            (2, 5),
            (2, 6),
            (3, 4),
            (3, 6),
            (4, 5),
            (5, 6)]
        """
        if not isinstance(adjacency,list):
            adjacency=[adjacency]
        pos_dict = {}
        for i in range(n):
            x = float(cos((pi/2) + ((2*pi)/n)*i))
            y = float(sin((pi/2) + ((2*pi)/n)*i))
            pos_dict[i] = (x,y)
        G=graph.Graph(n, name="Circulant graph ("+str(adjacency)+")")
        G._pos=pos_dict
        for v in G:
            G.add_edges([(v,(v+j)%n) for j in adjacency])
            G.add_edges([(v,(v-j)%n) for j in adjacency])
        return G



    def CompleteGraph(self, n):
        """
        Returns a complete graph on n nodes.

        A Complete Graph is a graph in which all nodes are connected to all
        other nodes.

        This constructor is dependent on vertices numbered 0 through n-1 in
        NetworkX complete_graph()

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, each complete
        graph will be displayed with the first (0) node at the top, with
        the rest following in a counterclockwise manner.

        In the complete graph, there is a big difference visually in using
        the spring-layout algorithm vs. the position dictionary used in
        this constructor. The position dictionary flattens the graph,
        making it clear which nodes an edge is connected to. But the
        complete graph offers a good example of how the spring-layout
        works. The edges push outward (everything is connected), causing
        the graph to appear as a 3-dimensional pointy ball. (See examples
        below).

        EXAMPLES: We view many Complete graphs with a Sage Graphics Array,
        first with this constructor (i.e., the position dictionary
        filled)::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.CompleteGraph(i+3)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time

        We compare to plotting with the spring-layout algorithm::

            sage: import networkx
            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    spr = networkx.complete_graph(i+3)
            ...    k = Graph(spr)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time

        Compare the constructors (results will vary)

        ::

            sage: import networkx
            sage: t = cputime()
            sage: n = networkx.complete_graph(389); spring389 = Graph(n)
            sage: cputime(t)           # random
            0.59203700000000126
            sage: t = cputime()
            sage: posdict389 = graphs.CompleteGraph(389)
            sage: cputime(t)           # random
            0.6680419999999998

        We compare plotting::

            sage: import networkx
            sage: n = networkx.complete_graph(23)
            sage: spring23 = Graph(n)
            sage: posdict23 = graphs.CompleteGraph(23)
            sage: spring23.show() # long time
            sage: posdict23.show() # long time
        """
        pos_dict = {}
        for i in range(n):
            x = float(cos((pi/2) + ((2*pi)/n)*i))
            y = float(sin((pi/2) + ((2*pi)/n)*i))
            pos_dict[i] = (x,y)
        import networkx
        G = networkx.complete_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Complete graph")

    def CompleteBipartiteGraph(self, n1, n2):
        """
        Returns a Complete Bipartite Graph sized n1+n2, with each of the
        nodes [0,(n1-1)] connected to each of the nodes [n1,(n2-1)] and
        vice versa.

        A Complete Bipartite Graph is a graph with its vertices partitioned
        into two groups, V1 and V2. Each v in V1 is connected to every v in
        V2, and vice versa.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, each complete
        bipartite graph will be displayed with the first n1 nodes on the
        top row (at y=1) from left to right. The remaining n2 nodes appear
        at y=0, also from left to right. The shorter row (partition with
        fewer nodes) is stretched to the same length as the longer row,
        unless the shorter row has 1 node; in which case it is centered.
        The x values in the plot are in domain [0,maxn1,n2].

        In the Complete Bipartite graph, there is a visual difference in
        using the spring-layout algorithm vs. the position dictionary used
        in this constructor. The position dictionary flattens the graph and
        separates the partitioned nodes, making it clear which nodes an
        edge is connected to. The Complete Bipartite graph plotted with the
        spring-layout algorithm tends to center the nodes in n1 (see
        spring_med in examples below), thus overlapping its nodes and
        edges, making it typically hard to decipher.

        Filling the position dictionary in advance adds O(n) to the
        constructor. Feel free to race the constructors below in the
        examples section. The much larger difference is the time added by
        the spring-layout algorithm when plotting. (Also shown in the
        example below). The spring model is typically described as
        `O(n^3)`, as appears to be the case in the NetworkX source
        code.

        EXAMPLES: Two ways of constructing the complete bipartite graph,
        using different layout algorithms::

            sage: import networkx
            sage: n = networkx.complete_bipartite_graph(389,157); spring_big = Graph(n)   # long time
            sage: posdict_big = graphs.CompleteBipartiteGraph(389,157)                    # long time

        Compare the plotting::

            sage: n = networkx.complete_bipartite_graph(11,17)
            sage: spring_med = Graph(n)
            sage: posdict_med = graphs.CompleteBipartiteGraph(11,17)

        Notice here how the spring-layout tends to center the nodes of n1

        ::

            sage: spring_med.show() # long time
            sage: posdict_med.show() # long time

        View many complete bipartite graphs with a Sage Graphics Array,
        with this constructor (i.e., the position dictionary filled)::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.CompleteBipartiteGraph(i+1,4)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time

        We compare to plotting with the spring-layout algorithm::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    spr = networkx.complete_bipartite_graph(i+1,4)
            ...    k = Graph(spr)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time
        """
        pos_dict = {}
        c1 = 1 # scaling factor for top row
        c2 = 1 # scaling factor for bottom row
        c3 = 0 # pad to center if top row has 1 node
        c4 = 0 # pad to center if bottom row has 1 node
        if n1 > n2:
            if n2 == 1:
                c4 = (n1-1)/2
            else:
                c2 = ((n1-1)/(n2-1))
        elif n2 > n1:
            if n1 == 1:
                c3 = (n2-1)/2
            else:
                c1 = ((n2-1)/(n1-1))
        for i in range(n1):
            x = c1*i + c3
            y = 1
            pos_dict[i] = (x,y)
        for i in range(n1+n2)[n1:]:
            x = c2*(i-n1) + c4
            y = 0
            pos_dict[i] = (x,y)
        import networkx
        import sage.graphs.bipartite_graph as bipartite_graph
        G = networkx.complete_bipartite_graph(n1,n2)
        return bipartite_graph.BipartiteGraph(G, pos=pos_dict, name="Complete bipartite graph")

    def CompleteMultipartiteGraph(self, l):
        r"""
        Returns a complete multipartite graph.

        INPUT:

        - ``l`` -- a list of integers : the respective sizes
          of the components.

        EXAMPLE:

        A complete tripartite graph with sets of sizes
        `5, 6, 8`::

            sage: g = graphs.CompleteMultipartiteGraph([5, 6, 8]); g
            Multipartite Graph with set sizes [5, 6, 8]: Graph on 19 vertices

        It clearly has a chromatic number of 3::

            sage: g.chromatic_number()
            3
        """

        from sage.graphs.graph import Graph
        g = Graph()
        for i in l:
            g = g + self.CompleteGraph(i)

        g = g.complement()
        g.name("Multipartite Graph with set sizes "+str(l))

        return g

    def CubeGraph(self, n):
        r"""
        Returns the hypercube in `n` dimensions.

        The hypercube in `n` dimension is build upon the binary
        strings on `n` bits, two of them being adjacent if
        they differ in exactly one bit. Hence, the distance
        between two vertices in the hypercube is the Hamming
        distance.

        EXAMPLES:

        The distance between `0100110` and `1011010` is
        `5`, as expected ::

            sage: g = graphs.CubeGraph(7)
            sage: g.distance('0100110','1011010')
            5

        Plot several `n`-cubes in a Sage Graphics Array ::

            sage: g = []
            sage: j = []
            sage: for i in range(6):
            ...    k = graphs.CubeGraph(i+1)
            ...    g.append(k)
            ...
            sage: for i in range(2):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show(figsize=[6,4]) # long time

        Use the plot options to display larger `n`-cubes

        ::

            sage: g = graphs.CubeGraph(9)
            sage: g.show(figsize=[12,12],vertex_labels=False, vertex_size=20) # long time

        AUTHORS:

        - Robert Miller
        """
        theta = float(pi/n)

        d = {'':[]}
        dn={}
        p = {'':(float(0),float(0))}
        pn={}

        # construct recursively the adjacency dict and the positions
        for i in xrange(n):
            ci = float(cos(i*theta))
            si = float(sin(i*theta))
            for v,e in d.iteritems():
                v0 = v+'0'
                v1 = v+'1'
                l0 = [v1]
                l1 = [v0]
                for m in e:
                    l0.append(m+'0')
                    l1.append(m+'1')
                dn[v0] = l0
                dn[v1] = l1
                x,y = p[v]
                pn[v0] = (x, y)
                pn[v1] = (x+ci, y+si)
            d,dn = dn,{}
            p,pn = pn,{}

        # construct the graph
        r = graph.Graph(name="%d-Cube"%n)
        r.add_vertices(d.keys())
        for u,L in d.iteritems():
            for v in L:
                r.add_edge(u,v)
        r.set_pos(p)

        return r

    def FuzzyBallGraph(self, partition, q):
        r"""
        Construct a Fuzzy Ball graph with the integer partition
        ``partition`` and ``q`` extra vertices.

        Let `q` be an integer and let `m_1,m_2,...,m_k` be a set of positive
        integers.  Let `n=q+m_1+...+m_k`.  The Fuzzy Ball graph with partition
        `m_1,m_2,...,m_k` and `q` extra vertices is the graph constructed from the
        graph `G=K_n` by attaching, for each `i=1,2,...,k`, a new vertex `a_i` to
        `m_i` distinct vertices of `G`.

        For given positive integers `k` and `m` and nonnegative
        integer `q`, the set of graphs ``FuzzyBallGraph(p, q)`` for
        all partitions `p` of `m` with `k` parts are cospectral with
        respect to the normalized Laplacian.

        EXAMPLES::

            sage: graphs.FuzzyBallGraph([3,1],2).adjacency_matrix()
            [0 1 1 1 1 1 1 0]
            [1 0 1 1 1 1 1 0]
            [1 1 0 1 1 1 1 0]
            [1 1 1 0 1 1 0 1]
            [1 1 1 1 0 1 0 0]
            [1 1 1 1 1 0 0 0]
            [1 1 1 0 0 0 0 0]
            [0 0 0 1 0 0 0 0]


        Pick positive integers `m` and `k` and a nonnegative integer `q`.
        All the FuzzyBallGraphs constructed from partitions of `m` with
        `k` parts should be cospectral with respect to the normalized
        Laplacian::

            sage: m=4; q=2; k=2
            sage: g_list=[graphs.FuzzyBallGraph(p,q) for p in Partitions(m, length=k)]
            sage: charpolys=set([g.laplacian_matrix(normalized=True).charpoly() for g in g_list])
            sage: charpolys
            set([x^8 - 8*x^7 + 4079/150*x^6 - 68689/1350*x^5 + 610783/10800*x^4 - 120877/3240*x^3 + 1351/100*x^2 - 931/450*x])
        """
        if len(partition)<1:
            raise ValueError, "partition must be a nonempty list of positive integers"
        n=q+sum(partition)
        g=graphs.CompleteGraph(n)
        curr_vertex=0
        for e,p in enumerate(partition):
            g.add_edges([(curr_vertex+i, 'a{0}'.format(e+1)) for i in range(p)])
            curr_vertex+=p
        return g



    def FibonacciTree(self, n):
        r"""
        Returns the graph of the Fibonacci Tree `F_{i}` of order `n`.
        `F_{i}` is recursively defined as the a tree with a root vertex
        and two attached child trees `F_{i-1}` and `F_{i-2}`, where
        `F_{1}` is just one vertex and `F_{0}` is empty.

        INPUT:

        - ``n`` - the recursion depth of the Fibonacci Tree

        EXAMPLES:

            sage: g = graphs.FibonacciTree(3)
            sage: g.is_tree()
            True
            sage: l1 = [ len(graphs.FibonacciTree(_)) + 1 for _ in range(6) ]
            sage: l2 = list(fibonacci_sequence(2,8))
            sage: l1 == l2
            True

        AUTHORS:

        - Harald Schilly and Yann Laigle-Chapuy (2010-03-25)
        """
        T = graph.Graph(name="Fibonacci-Tree-%d"%n)
        if n == 1: T.add_vertex(0)
        if n < 2: return T

        from sage.combinat.combinat import fibonacci_sequence
        F = list(fibonacci_sequence(n + 2))
        s = 1.618 ** (n / 1.618 - 1.618)
        pos = {}

        def fib(level, node, y):
            pos[node] = (node, y)
            if level < 2: return
            level -= 1
            y -= s
            diff = F[level]
            T.add_edge(node, node - diff)
            if level == 1: # only one child
                pos[node - diff] = (node, y)
                return
            T.add_edge(node, node + diff)
            fib(level, node - diff, y)
            fib(level - 1, node + diff, y)

        T.add_vertices(xrange(sum(F[:-1])))
        fib(n, F[n + 1] - 1, 0)
        T.set_pos(pos)

        return T

    def GeneralizedPetersenGraph(self, n,k):
        r"""
        Returns a generalized Petersen graph with `2n` nodes. The variables
        `n`, `k` are integers such that `n>2` and `0<k\leq\lfloor(n-1)`/`2\rfloor`

        For `k=1` the result is a graph isomorphic to the circular ladder graph
        with the same `n`. The regular Petersen Graph has `n=5` and `k=2`.
        Other named graphs that can be described using this notation include
        the Desargues graph and the Moebius-Kantor graph.

        INPUT:

        - ``n`` - the number of nodes is `2*n`.

        - ``k`` - integer `0<k\leq\lfloor(n-1)`/`2\rfloor`. Decides
          how inner vertices are connected.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the generalized
        Petersen graphs are displayed as an inner and outer cycle pair, with
        the first n nodes drawn on the outer circle. The first (0) node is
        drawn at the top of the outer-circle, moving counterclockwise after that.
        The inner circle is drawn with the (n)th node at the top, then
        counterclockwise as well.

        EXAMPLES: For `k=1` the resulting graph will be isomorphic to a circular
        ladder graph. ::

            sage: g = graphs.GeneralizedPetersenGraph(13,1)
            sage: g2 = graphs.CircularLadderGraph(13)
            sage: g.is_isomorphic(g2)
            True

        The Desargues graph::

            sage: g = graphs.GeneralizedPetersenGraph(10,3)
            sage: g.girth()
            6
            sage: g.is_bipartite()
            True

        AUTHORS:

        - Anders Jonsson (2009-10-15)
        """
        if (n < 3):
                raise ValueError("n must be larger than 2")
        if (k < 1 or k>((n-1)/2)):
                raise ValueError("k must be in 1<= k <=floor((n-1)/2)")
        pos_dict = {}
        G=graphs.EmptyGraph()
        for i in range(n):
            x = float(cos((pi/2) + ((2*pi)/n)*i))
            y = float(sin((pi/2) + ((2*pi)/n)*i))
            pos_dict[i] = (x,y)
        for i in range(n, 2*n):
            x = float(0.5*cos((pi/2) + ((2*pi)/n)*i))
            y = float(0.5*sin((pi/2) + ((2*pi)/n)*i))
            pos_dict[i] = (x,y)
        for i in range(n):
            G.add_edge(i, (i+1) % n)
            G.add_edge(i, i+n)
            G.add_edge(i+n, n + (i+k) % n)
        return graph.Graph(G, pos=pos_dict, name="Generalized Petersen graph (n="+str(n)+",k="+str(k)+")")

    def HyperStarGraph(self,n,k):
        r"""
        Returns the hyper-star graph HS(n,k).

        The vertices of the hyper-star graph are the set of binary strings
        of length n which contain k 1s. Two vertices, u and v, are adjacent
        only if u can be obtained from v by swapping the first bit with a
        different symbol in another position.

        INPUT:

        -  ``n``

        -  ``k``

        EXAMPLES::

            sage: g = graphs.HyperStarGraph(6,3)
            sage: g.plot() # long time

        REFERENCES:

        - Lee, Hyeong-Ok, Jong-Seok Kim, Eunseuk Oh, and Hyeong-Seok Lim.
          "Hyper-Star Graph: A New Interconnection Network Improving the
          Network Cost of the Hypercube." In Proceedings of the First EurAsian
          Conference on Information and Communication Technology, 858-865.
          Springer-Verlag, 2002.

        AUTHORS:

        - Michael Yurko (2009-09-01)
        """
        from sage.combinat.combination import Combinations
        # dictionary associating the positions of the 1s to the corresponding
        # string: e.g. if n=6 and k=3, comb_to_str([0,1,4])=='110010'
        comb_to_str={}
        for c in Combinations(n,k):
            L = ['0']*n
            for i in c:
                L[i]='1'
            comb_to_str[tuple(c)] = ''.join(L)

        g=graph.Graph(name="HS(%d,%d)"%(n,k))
        g.add_vertices(comb_to_str.values())

        for c in Combinations(range(1,n),k): # 0 is not in c
            L = []
            u = comb_to_str[tuple(c)]
            # switch 0 with the 1s
            for i in xrange(len(c)):
                v = tuple([0]+c[:i]+c[i+1:])
                g.add_edge( u , comb_to_str[v] )

        return g

    def NKStarGraph(self,n,k):
        r"""
        Returns the (n,k)-star graph.

        The vertices of the (n,k)-star graph are the set of all arrangements of
        n symbols into labels of length k. There are two adjacency rules for
        the (n,k)-star graph. First, two vertices are adjacent if one can be
        obtained from the other by swapping the first symbol with another
        symbol. Second, two vertices are adjacent if one can be obtained from
        the other by swapping the first symbol with an external symbol (a
        symbol not used in the original label).

        INPUT:

        -  ``n``

        -  ``k``

        EXAMPLES::

            sage: g = graphs.NKStarGraph(4,2)
            sage: g.plot() # long time

        REFERENCES:

        - Wei-Kuo, Chiang, and Chen Rong-Jaye. "The (n, k)-star graph: A
          generalized star graph." Information Processing Letters 56,
          no. 5 (December 8, 1995): 259-264.

        AUTHORS:

        - Michael Yurko (2009-09-01)
        """
        from sage.combinat.permutation import Arrangements
        #set from which to permute
        set = [str(i) for i in xrange(1,n+1)]
        #create dict
        d = {}
        for v in Arrangements(set,k):
            tmp_dict = {}
            #add edges of dimension i
            for i in xrange(1,k):
                #swap 0th and ith element
                v[0], v[i] = v[i], v[0]
                #convert to str and add to list
                vert = "".join(v)
                tmp_dict[vert] = None
                #swap back
                v[0], v[i] = v[i], v[0]
            #add other edges
            tmp_bit = v[0]
            for i in set:
                #check if external
                if not (i in v):
                    v[0] = i
                    #add edge
                    vert = "".join(v)
                    tmp_dict[vert] = None
                v[0] = tmp_bit
            d["".join(v)] = tmp_dict
        return graph.Graph(d, name="(%d,%d)-star"%(n,k))

    def NStarGraph(self,n):
        r"""
        Returns the n-star graph.

        The vertices of the n-star graph are the set of permutations on n
        symbols. There is an edge between two vertices if their labels differ
        only in the first and one other position.

        INPUT:

        -  ``n``

        EXAMPLES::

            sage: g = graphs.NStarGraph(4)
            sage: g.plot() # long time

        REFERENCES:

        - S.B. Akers, D. Horel and B. Krishnamurthy, The star graph: An
          attractive alternative to the previous n-cube. In: Proc. Internat.
          Conf. on Parallel Processing (1987), pp. 393--400.

        AUTHORS:

        - Michael Yurko (2009-09-01)
        """
        from sage.combinat.permutation import Permutations
        #set from which to permute
        set = [str(i) for i in xrange(1,n+1)]
        #create dictionary of lists
        #vertices are adjacent if the first element
        #is swapped with the ith element
        d = {}
        for v in Permutations(set):
            tmp_dict = {}
            for i in xrange(1,n):
                if v[0] != v[i]:
                    #swap 0th and ith element
                    v[0], v[i] = v[i], v[0]
                    #convert to str and add to list
                    vert = "".join(v)
                    tmp_dict[vert] = None
                    #swap back
                    v[0], v[i] = v[i], v[0]
            d["".join(v)] = tmp_dict
        return graph.Graph(d, name = "%d-star"%n)

    def BubbleSortGraph(self, n):
        r"""
        Returns the bubble sort graph `B(n)`.

        The vertices of the bubble sort graph are the set of permutations on
        `n` symbols. Two vertices are adjacent if one can be obtained from the
        other by swapping the labels in the `i`-th and `(i+1)`-th positions for
        `1 \leq i \leq n-1`. In total, `B(n)` has order `n!`. Thus, the order
        of `B(n)` increases according to `f(n) = n!`.

        INPUT:

        - ``n`` -- positive integer. The number of symbols to permute.

        OUTPUT:

        The bubble sort graph `B(n)` on `n` symbols. If `n < 1`, a
        ``ValueError`` is returned.

        EXAMPLES::

            sage: g = graphs.BubbleSortGraph(4); g
            Bubble sort: Graph on 24 vertices
            sage: g.plot() # long time

        The bubble sort graph on `n = 1` symbol is the trivial graph `K_1`::

            sage: graphs.BubbleSortGraph(1)
            Bubble sort: Graph on 1 vertex

        If `n \geq 1`, then the order of `B(n)` is `n!`::

            sage: n = randint(1, 8)
            sage: g = graphs.BubbleSortGraph(n)
            sage: g.order() == factorial(n)
            True

        TESTS:

        Input ``n`` must be positive::

            sage: graphs.BubbleSortGraph(0)
            Traceback (most recent call last):
            ...
            ValueError: Invalid number of symbols to permute, n should be >= 1
            sage: graphs.BubbleSortGraph(randint(-10^6, 0))
            Traceback (most recent call last):
            ...
            ValueError: Invalid number of symbols to permute, n should be >= 1

        AUTHORS:

        - Michael Yurko (2009-09-01)
        """
        # sanity checks
        if n < 1:
            raise ValueError(
                "Invalid number of symbols to permute, n should be >= 1")
        if n == 1:
            return graph.Graph(self.CompleteGraph(n), name="Bubble sort")
        from sage.combinat.permutation import Permutations
        #create set from which to permute
        label_set = [str(i) for i in xrange(1, n + 1)]
        d = {}
        #iterate through all vertices
        for v in Permutations(label_set):
            tmp_dict = {}
            #add all adjacencies
            for i in xrange(n - 1):
                #swap entries
                v[i], v[i + 1] = v[i + 1], v[i]
                #add new vertex
                new_vert = ''.join(v)
                tmp_dict[new_vert] = None
                #swap back
                v[i], v[i + 1] = v[i + 1], v[i]
            #add adjacency dict
            d[''.join(v)] = tmp_dict
        return graph.Graph(d, name="Bubble sort")

    def BalancedTree(self, r, h):
        r"""
        Returns the perfectly balanced tree of height `h \geq 1`,
        whose root has degree `r \geq 2`.

        The number of vertices of this graph is
        `1 + r + r^2 + \cdots + r^h`, that is,
        `\frac{r^{h+1} - 1}{r - 1}`. The number of edges is one
        less than the number of vertices.

        INPUT:

        - ``r`` -- positive integer `\geq 2`. The degree of the root node.

        - ``h`` -- positive integer `\geq 1`. The height of the balanced tree.

        OUTPUT:

        The perfectly balanced tree of height `h \geq 1` and whose root has
        degree `r \geq 2`. A ``NetworkXError`` is returned if `r < 2` or
        `h < 1`.

        ALGORITHM:

        Uses `NetworkX <http://networkx.lanl.gov>`_.

        EXAMPLES:

        A balanced tree whose root node has degree `r = 2`, and of height
        `h = 1`, has order 3 and size 2::

            sage: G = graphs.BalancedTree(2, 1); G
            Balanced tree: Graph on 3 vertices
            sage: G.order(); G.size()
            3
            2
            sage: r = 2; h = 1
            sage: v = 1 + r
            sage: v; v - 1
            3
            2

        Plot a balanced tree of height 5, whose root node has degree `r = 3`::

            sage: G = graphs.BalancedTree(3, 5)
            sage: G.show()   # long time

        A tree is bipartite. If its vertex set is finite, then it is planar. ::

            sage: r = randint(2, 5); h = randint(1, 7)
            sage: T = graphs.BalancedTree(r, h)
            sage: T.is_bipartite()
            True
            sage: T.is_planar()
            True
            sage: v = (r^(h + 1) - 1) / (r - 1)
            sage: T.order() == v
            True
            sage: T.size() == v - 1
            True

        TESTS:

        We only consider balanced trees whose root node has degree `r \geq 2`::

            sage: graphs.BalancedTree(1, randint(1, 10^6))
            Traceback (most recent call last):
            ...
            NetworkXError: Invalid graph description, r should be >=2
            sage: graphs.BalancedTree(randint(-10^6, 1), randint(1, 10^6))
            Traceback (most recent call last):
            ...
            NetworkXError: Invalid graph description, r should be >=2

        The tree must have height `h \geq 1`::

            sage: graphs.BalancedTree(randint(2, 10^6), 0)
            Traceback (most recent call last):
            ...
            NetworkXError: Invalid graph description, h should be >=1
            sage: graphs.BalancedTree(randint(2, 10^6), randint(-10^6, 0))
            Traceback (most recent call last):
            ...
            NetworkXError: Invalid graph description, h should be >=1
            sage: graphs.BalancedTree(randint(-10^6, 1), randint(-10^6, 0))
            Traceback (most recent call last):
            ...
            NetworkXError: Invalid graph description, r should be >=2
        """
        import networkx
        return graph.Graph(networkx.balanced_tree(r, h), name="Balanced tree")

    def LCFGraph(self, n, shift_list, repeats):
        """
        Returns the cubic graph specified in LCF notation.

        LCF (Lederberg-Coxeter-Fruchte) notation is a concise way of
        describing cubic Hamiltonian graphs. The way a graph is constructed
        is as follows. Since there is a Hamiltonian cycle, we first create
        a cycle on n nodes. The variable shift_list = [s_0, s_1, ...,
        s_k-1] describes edges to be created by the following scheme: for
        each i, connect vertex i to vertex (i + s_i). Then, repeats
        specifies the number of times to repeat this process, where on the
        jth repeat we connect vertex (i + j\*len(shift_list)) to vertex (
        i + j\*len(shift_list) + s_i).

        INPUT:


        -  ``n`` - the number of nodes.

        -  ``shift_list`` - a list of integer shifts mod n.

        -  ``repeats`` - the number of times to repeat the
           process.


        EXAMPLES::

            sage: G = graphs.LCFGraph(4, [2,-2], 2)
            sage: G.is_isomorphic(graphs.TetrahedralGraph())
            True

        ::

            sage: G = graphs.LCFGraph(20, [10,7,4,-4,-7,10,-4,7,-7,4], 2)
            sage: G.is_isomorphic(graphs.DodecahedralGraph())
            True

        ::

            sage: G = graphs.LCFGraph(14, [5,-5], 7)
            sage: G.is_isomorphic(graphs.HeawoodGraph())
            True

        The largest cubic nonplanar graph of diameter three::

            sage: G = graphs.LCFGraph(20, [-10,-7,-5,4,7,-10,-7,-4,5,7,-10,-7,6,-5,7,-10,-7,5,-6,7], 1)
            sage: G.degree()
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
            sage: G.diameter()
            3
            sage: G.show()  # long time

        PLOTTING: LCF Graphs are plotted as an n-cycle with edges in the
        middle, as described above.

        REFERENCES:

        - [1] Frucht, R. "A Canonical Representation of Trivalent
          Hamiltonian Graphs." J. Graph Th. 1, 45-60, 1976.

        - [2] Grunbaum, B.  Convex Polytope es. New York: Wiley,
          pp. 362-364, 1967.

        - [3] Lederberg, J. 'DENDRAL-64: A System for Computer
          Construction, Enumeration and Notation of Organic Molecules
          as Tree Structures and Cyclic Graphs. Part II. Topology of
          Cyclic Graphs.' Interim Report to the National Aeronautics
          and Space Administration. Grant NsG 81-60. December 15,
          1965.  http://profiles.nlm.nih.gov/BB/A/B/I/U/_/bbabiu.pdf.
        """
        import networkx
        pos_dict = {}
        for i in range(n):
            x = float(cos(pi/2 + ((2*pi)/n)*i))
            y = float(sin(pi/2 + ((2*pi)/n)*i))
            pos_dict[i] = [x,y]
        return graph.Graph(networkx.LCF_graph(n, shift_list, repeats),\
                           pos=pos_dict, name="LCF Graph")

################################################################################
#   Pseudofractal Graphs
################################################################################

    def DorogovtsevGoltsevMendesGraph(self, n):
        """
        Construct the n-th generation of the Dorogovtsev-Goltsev-Mendes
        graph.

        EXAMPLE::

            sage: G = graphs.DorogovtsevGoltsevMendesGraph(8)
            sage: G.size()
            6561

        REFERENCE:

        - [1] Dorogovtsev, S. N., Goltsev, A. V., and Mendes, J.
          F. F., Pseudofractal scale-free web, Phys. Rev. E 066122
          (2002).
        """
        import networkx
        return graph.Graph(networkx.dorogovtsev_goltsev_mendes_graph(n),\
               name="Dorogovtsev-Goltsev-Mendes Graph, %d-th generation"%n)

################################################################################
#   Random Graphs
################################################################################

    def RandomGNP(self, n, p, seed=None, fast=True):
        r"""
        Returns a Random graph on `n` nodes. Each edge is inserted
        independently with probability `p`.

        IMPLEMENTATION: This function calls the NetworkX function
        ``fast_gnp_random_graph``, unless fast==False, then
        ``gnp_random_graph``.

        REFERENCES:

        - [1] P. Erdos and A. Renyi, On Random Graphs, Publ.
          Math. 6, 290 (1959).

        - [2] E. N. Gilbert, Random Graphs, Ann. Math.
          Stat., 30, 1141 (1959).

        PLOTTING: When plotting, this graph will use the default
        spring-layout algorithm, unless a position dictionary is
        specified.

        EXAMPLES: We show the edge list of a random graph on 6 nodes with
        probability `p = .4`::

            sage: graphs.RandomGNP(6, .4).edges(labels=False)
            [(0, 1), (0, 3), (0, 4), (0, 5), (1, 2), (1, 3), (1, 4), (1, 5)]

        We plot a random graph on 12 nodes with probability
        `p = .71`::

            sage: gnp = graphs.RandomGNP(12,.71)
            sage: gnp.show() # long time

        We view many random graphs using a graphics array::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.RandomGNP(i+3,.43)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show() # long time
            sage: graphs.RandomGNP(4,1)
            Complete graph: Graph on 4 vertices

        TIMINGS: The following timings compare the speed with fast==False
        and fast==True for sparse and dense graphs. (It's no different?)

        ::

            sage: t=cputime(); regular_sparse = graphs.RandomGNP(389,.22)
            sage: cputime(t)     # slightly random
            0.2240130000000029

        ::

            sage: t=cputime(); fast_sparse =  graphs.RandomGNP(389,.22,fast=True)
            sage: cputime(t)     # slightly random
            0.22401400000000038

        ::

            sage: t=cputime(); regular_dense = graphs.RandomGNP(389,.88)    # long time
            sage: cputime(t)     # slightly random, long time
            0.87205499999999958

        ::

            sage: t=cputime(); fast_dense = graphs.RandomGNP(389,.88,fast=True)    # long time
            sage: cputime(t)     # slightly random, long time
            0.90005700000000033
        """
        if seed is None:
            seed = current_randstate().long_seed()
        if p == 1:
            return graphs.CompleteGraph(n)
        import networkx
        if fast:
            G = networkx.fast_gnp_random_graph(n, p, seed=seed)
        else:
            G = networkx.gnp_random_graph(n, p, seed=seed)
        return graph.Graph(G)

    def RandomBarabasiAlbert(self, n, m, seed=None):
        u"""
        Return a random graph created using the Barabasi-Albert preferential
        attachment model.

        A graph with m vertices and no edges is initialized, and a graph of n
        vertices is grown by attaching new vertices each with m edges that are
        attached to existing vertices, preferentially with high degree.

        INPUT:

        - ``n`` - number of vertices in the graph

        - ``m`` - number of edges to attach from each new node

        - ``seed`` - for random number generator

        EXAMPLES:

        We show the edge list of a random graph on 6 nodes with m = 2.

        ::

            sage: graphs.RandomBarabasiAlbert(6,2).edges(labels=False)
            [(0, 2), (0, 3), (0, 4), (1, 2), (2, 3), (2, 4), (2, 5), (3, 5)]

        We plot a random graph on 12 nodes with m = 3.

        ::

            sage: ba = graphs.RandomBarabasiAlbert(12,3)
            sage: ba.show()  # long time

        We view many random graphs using a graphics array::

            sage: g = []
            sage: j = []
            sage: for i in range(1,10):
            ...    k = graphs.RandomBarabasiAlbert(i+3, 3)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show()  # long time

        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return graph.Graph(networkx.barabasi_albert_graph(n,m,seed=seed))

    def RandomBipartite(self, n1,n2, p):
        r"""
        Returns a bipartite graph with `n1+n2` vertices
        such that any edge from `[n1]` to `[n2]` exists
        with probability `p`.

        INPUT:

            - ``n1,n2`` : Cardinalities of the two sets
            - ``p``   : Probability for an edge to exist


        EXAMPLE::

            sage: g=graphs.RandomBipartite(5,2,0.5)
            sage: g.vertices()
            [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (1, 0), (1, 1)]

        TESTS::

            sage: g=graphs.RandomBipartite(5,-3,0.5)
            Traceback (most recent call last):
            ...
            ValueError: n1 and n2 should be integers strictly greater than 0
            sage: g=graphs.RandomBipartite(5,3,1.5)
            Traceback (most recent call last):
            ...
            ValueError: Parameter p is a probability, and so should be a real value between 0 and 1

        """
        if not (p>=0 and p<=1):
            raise ValueError, "Parameter p is a probability, and so should be a real value between 0 and 1"
        if not (n1>0 and n2>0):
            raise ValueError, "n1 and n2 should be integers strictly greater than 0"

        from numpy.random import uniform
        import sage.graphs.bipartite_graph as bipartite_graph
        from sage.graphs.all import Graph

        g=Graph()

        S1=[(0,i) for i in range(n1)]
        S2=[(1,i) for i in range(n2)]
        g.add_vertices(S1)
        g.add_vertices(S2)
        [g.add_edge((0,v),(1,w)) for v in range(n1) for w in range(n2) if uniform()<=p]

        return bipartite_graph.BipartiteGraph(g,[S1,S2],name="Random bipartite graph of size "+str(n1) +"+"+str(n2)+" with edge probability "+str(p))

    def RandomGNM(self, n, m, dense=False, seed=None):
        """
        Returns a graph randomly picked out of all graphs on n vertices
        with m edges.

        INPUT:


        -  ``n`` - number of vertices.

        -  ``m`` - number of edges.

        -  ``dense`` - whether to use NetworkX's
           dense_gnm_random_graph or gnm_random_graph


        EXAMPLES: We show the edge list of a random graph on 5 nodes with
        10 edges.

        ::

            sage: graphs.RandomGNM(5, 10).edges(labels=False)
            [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]

        We plot a random graph on 12 nodes with m = 12.

        ::

            sage: gnm = graphs.RandomGNM(12, 12)
            sage: gnm.show()  # long time

        We view many random graphs using a graphics array::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.RandomGNM(i+3, i^2-i)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.show()  # long time
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        if dense:
            return graph.Graph(networkx.dense_gnm_random_graph(n, m, seed=seed))
        else:
            return graph.Graph(networkx.gnm_random_graph(n, m, seed=seed))

    def RandomNewmanWattsStrogatz(self, n, k, p, seed=None):
        """
        Returns a Newman-Watts-Strogatz small world random graph on n
        vertices.

        From the NetworkX documentation: First create a ring over n nodes.
        Then each node in the ring is connected with its k nearest
        neighbors. Then shortcuts are created by adding new edges as
        follows: for each edge u-v in the underlying "n-ring with k nearest
        neighbors"; with probability p add a new edge u-w with
        randomly-chosen existing node w. In contrast with
        watts_strogatz_graph(), no edges are removed.

        INPUT:


        -  ``n`` - number of vertices.

        -  ``k`` - each vertex is connected to its k nearest
           neighbors

        -  ``p`` - the probability of adding a new edge for
           each edge

        -  ``seed`` - for the random number generator


        EXAMPLE: We show the edge list of a random graph on 7 nodes with 2
        "nearest neighbors" and probability `p = 0.2`::

            sage: graphs.RandomNewmanWattsStrogatz(7, 2, 0.2).edges(labels=False)
            [(0, 1), (0, 2), (0, 3), (0, 6), (1, 2), (2, 3), (2, 4), (3, 4), (3, 6), (4, 5), (5, 6)]

        ::

            sage: G = graphs.RandomNewmanWattsStrogatz(12, 2, .3)
            sage: G.show()  # long time

        REFERENCE:

        - [1] Newman, M.E.J., Watts, D.J. and Strogatz, S.H.  Random
          graph models of social networks. Proc. Nat. Acad. Sci. USA
          99, 2566-2572.
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return graph.Graph(networkx.newman_watts_strogatz_graph(n, k, p, seed=seed))

    def RandomHolmeKim(self, n, m, p, seed=None):
        """
        Returns a random graph generated by the Holme and Kim algorithm for
        graphs with power law degree distribution and approximate average
        clustering.

        INPUT:


        -  ``n`` - number of vertices.

        -  ``m`` - number of random edges to add for each new
           node.

        -  ``p`` - probability of adding a triangle after
           adding a random edge.

        -  ``seed`` - for the random number generator.


        From the NetworkX documentation: The average clustering has a hard
        time getting above a certain cutoff that depends on m. This cutoff
        is often quite low. Note that the transitivity (fraction of
        triangles to possible triangles) seems to go down with network
        size. It is essentially the Barabasi-Albert growth model with an
        extra step that each random edge is followed by a chance of making
        an edge to one of its neighbors too (and thus a triangle). This
        algorithm improves on B-A in the sense that it enables a higher
        average clustering to be attained if desired. It seems possible to
        have a disconnected graph with this algorithm since the initial m
        nodes may not be all linked to a new node on the first iteration
        like the BA model.

        EXAMPLE: We show the edge list of a random graph on 8 nodes with 2
        random edges per node and a probability `p = 0.5` of
        forming triangles.

        ::

            sage: graphs.RandomHolmeKim(8, 2, 0.5).edges(labels=False)
            [(0, 2), (0, 4), (1, 2), (1, 3), (2, 3), (3, 4), (3, 5), (3, 6), (3, 7), (4, 5), (4, 6)]

        ::

            sage: G = graphs.RandomHolmeKim(12, 3, .3)
            sage: G.show()  # long time

        REFERENCE:

        - [1] Holme, P. and Kim, B.J. Growing scale-free networks with
          tunable clustering, Phys. Rev. E (2002). vol 65, no 2,
          026107.
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return graph.Graph(networkx.powerlaw_cluster_graph(n, m, p, seed=seed))

    def RandomInterval(self,n):
        """
        Returns a random interval graph.

        An interval graph is built from a list `(a_i,b_i)_{1\leq i \leq n}`
        of intervals : to each interval of the list is associated one
        vertex, two vertices being adjacent if the two corresponding
        intervals intersect.

        A random interval graph of order `n` is generated by picking
        random values for the `(a_i,b_j)`, each of the two coordinates
        being generated from the uniform distribution on the interval
        `[0,1]`.

        This definitions follows [boucheron2001]_.

        INPUT:

        - ``n`` (integer) -- the number of vertices in the random
          graph.

        EXAMPLE:

        As for any interval graph, the chromatic number is equal to
        the clique number ::

            sage: g = graphs.RandomInterval(8)
            sage: g.clique_number() == g.chromatic_number()
            True

        REFERENCE:

        .. [boucheron2001] Boucheron, S. and FERNANDEZ de la VEGA, W.,
           On the Independence Number of Random Interval Graphs,
           Combinatorics, Probability and Computing v10, issue 05,
           Pages 385--396,
           Cambridge Univ Press, 2001
        """

        from sage.misc.prandom import random

        intervals = [tuple(sorted((random(), random()))) for i in range(n)]

        return self.IntervalGraph(intervals)

    def IntervalGraph(self,intervals):
        r"""
        Returns the graph corresponding to the given intervals.

        An interval graph is built from a list `(a_i,b_i)_{1\leq i \leq n}`
        of intervals : to each interval of the list is associated one
        vertex, two vertices being adjacent if the two corresponding
        (closed) intervals intersect.

        INPUT:

        - ``intervals`` -- the list of pairs `(a_i,b_i)`
          defining the graph.

        NOTE:

        The intervals `(a_i,b_i)` must verify `a_i<b_i`.

        EXAMPLE:

        The following line creates the sequence of intervals
        `(i, i+2)` for i in `[0, ..., 8]`::

            sage: intervals = [(i,i+2) for i in range(9)]

        In the corresponding graph...::

            sage: g = graphs.IntervalGraph(intervals)
            sage: sorted(g.neighbors((3,5)))
            [(1, 3), (2, 4), (4, 6), (5, 7)]

        And the clique number is as expected ::

            sage: g.clique_number()
            3
        """

        intervals = sorted(intervals)


        edges = []
        while intervals:
            x = intervals.pop(0)
            for y in intervals:
                if y[0] <= x[1]:
                    edges.append((x,y))
                else:
                    break
        g = graph.Graph(vertices=intervals)
        g.add_edges(edges)
        return g

    def RandomLobster(self, n, p, q, seed=None):
        """
        Returns a random lobster.

        A lobster is a tree that reduces to a caterpillar when pruning all
        leaf vertices. A caterpillar is a tree that reduces to a path when
        pruning all leaf vertices (q=0).

        INPUT:


        -  ``n`` - expected number of vertices in the backbone

        -  ``p`` - probability of adding an edge to the
           backbone

        -  ``q`` - probability of adding an edge (claw) to the
           arms

        -  ``seed`` - for the random number generator


        EXAMPLE: We show the edge list of a random graph with 3 backbone
        nodes and probabilities `p = 0.7` and `q = 0.3`::

            sage: graphs.RandomLobster(3, 0.7, 0.3).edges(labels=False)
            [(0, 1), (1, 2)]

        ::

            sage: G = graphs.RandomLobster(9, .6, .3)
            sage: G.show()  # long time
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return graph.Graph(networkx.random_lobster(n, p, q, seed=seed))

    def RandomTreePowerlaw(self, n, gamma=3, tries=100, seed=None):
        """
        Returns a tree with a power law degree distribution. Returns False
        on failure.

        From the NetworkX documentation: A trial power law degree sequence
        is chosen and then elements are swapped with new elements from a
        power law distribution until the sequence makes a tree (size = order
        - 1).

        INPUT:


        -  ``n`` - number of vertices

        -  ``gamma`` - exponent of power law

        -  ``tries`` - number of attempts to adjust sequence to
           make a tree

        -  ``seed`` - for the random number generator


        EXAMPLE: We show the edge list of a random graph with 10 nodes and
        a power law exponent of 2.

        ::

            sage: graphs.RandomTreePowerlaw(10, 2).edges(labels=False)
            [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (6, 8), (6, 9)]

        ::

            sage: G = graphs.RandomTreePowerlaw(15, 2)
            sage: if G:
            ...    G.show()  # random output, long time
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        try:
            return graph.Graph(networkx.random_powerlaw_tree(n, gamma, seed=seed, tries=tries))
        except:
            return False

    def RandomRegular(self, d, n, seed=None):
        """
        Returns a random d-regular graph on n vertices, or returns False on
        failure.

        Since every edge is incident to two vertices, n\*d must be even.

        INPUT:


        -  ``n`` - number of vertices

        -  ``d`` - degree

        -  ``seed`` - for the random number generator


        EXAMPLE: We show the edge list of a random graph with 8 nodes each
        of degree 3.

        ::

            sage: graphs.RandomRegular(3, 8).edges(labels=False)
            [(0, 1), (0, 4), (0, 7), (1, 5), (1, 7), (2, 3), (2, 5), (2, 6), (3, 4), (3, 6), (4, 5), (6, 7)]

        ::

            sage: G = graphs.RandomRegular(3, 20)
            sage: if G:
            ...    G.show()  # random output, long time

        REFERENCES:

        - [1] Kim, Jeong Han and Vu, Van H. Generating random regular
          graphs. Proc. 35th ACM Symp. on Thy. of Comp. 2003, pp
          213-222. ACM Press, San Diego, CA, USA.
          http://doi.acm.org/10.1145/780542.780576

        - [2] Steger, A. and Wormald, N. Generating random regular
          graphs quickly. Prob. and Comp. 8 (1999), pp 377-396.
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        try:
            N = networkx.random_regular_graph(d, n, seed=seed)
            if N is False: return False
            return graph.Graph(N, sparse=True)
        except:
            return False

    def RandomShell(self, constructor, seed=None):
        """
        Returns a random shell graph for the constructor given.

        INPUT:


        -  ``constructor`` - a list of 3-tuples (n,m,d), each
           representing a shell

        -  ``n`` - the number of vertices in the shell

        -  ``m`` - the number of edges in the shell

        -  ``d`` - the ratio of inter (next) shell edges to
           intra shell edges

        -  ``seed`` - for the random number generator


        EXAMPLE::

            sage: G = graphs.RandomShell([(10,20,0.8),(20,40,0.8)])
            sage: G.edges(labels=False)
            [(0, 3), (0, 7), (0, 8), (1, 2), (1, 5), (1, 8), (1, 9), (3, 6), (3, 11), (4, 6), (4, 7), (4, 8), (4, 21), (5, 8), (5, 9), (6, 9), (6, 10), (7, 8), (7, 9), (8, 18), (10, 11), (10, 13), (10, 19), (10, 22), (10, 26), (11, 18), (11, 26), (11, 28), (12, 13), (12, 14), (12, 28), (12, 29), (13, 16), (13, 21), (13, 29), (14, 18), (16, 20), (17, 18), (17, 26), (17, 28), (18, 19), (18, 22), (18, 27), (18, 28), (19, 23), (19, 25), (19, 28), (20, 22), (24, 26), (24, 27), (25, 27), (25, 29)]
            sage: G.show()  # long time
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return graph.Graph(networkx.random_shell_graph(constructor, seed=seed))

    def WorldMap(self):
        """
        Returns the Graph of all the countries, in which two countries are adjacent
        in the graph if they have a common boundary.

        This graph has been built from the data available
        in The CIA World Factbook [CIA]_ (2009-08-21).

        The returned graph ``G`` has a member ``G.gps_coordinates``
        equal to a dictionary containing the GPS coordinates
        of each country's capital city.

        EXAMPLE::

            sage: g=graphs.WorldMap()
            sage: g.has_edge("France","Italy")
            True
            sage: g.gps_coordinates["Bolivia"]
            [[17, 'S'], [65, 'W']]
            sage: sorted(g.connected_component_containing_vertex('Ireland'))
            ['Ireland', 'United Kingdom']

        REFERENCE:

        .. [CIA] CIA Factbook 09 https://www.cia.gov/library/publications/the-world-factbook/
        """
        from sage.graphs.all import Graph
        edges = [
            ('Afghanistan', 'China', None), ('Afghanistan', 'Iran', None),
            ('Afghanistan', 'Uzbekistan', None), ('Albania', 'Greece', None),
            ('Albania', 'Kosovo', None), ('Albania', 'Macedonia', None),
            ('Albania', 'Montenegro', None), ('Algeria', 'Morocco', None),
            ('Algeria', 'Tunisia', None), ('Andorra', 'Spain', None),
            ('Angola', 'Democratic Republic of the Congo', None), ('Angola', 'Namibia', None),
            ('Angola', 'Zambia', None), ('Argentina', 'Bolivia', None),
            ('Argentina', 'Brazil', None), ('Argentina', 'Chile', None),
            ('Argentina', 'Paraguay', None), ('Argentina', 'Uruguay', None),
            ('Armenia', 'Georgia', None), ('Armenia', 'Iran', None),
            ('Austria', 'Germany', None), ('Azerbaijan', 'Armenia', None),
            ('Azerbaijan', 'Georgia', None), ('Azerbaijan', 'Iran', None),
            ('Azerbaijan', 'Russia', None), ('Azerbaijan', 'Turkey', None),
            ('Bangladesh', 'Burma', None), ('Belgium', 'Germany', None),
            ('Belgium', 'Netherlands', None), ('Belize', 'Mexico', None),
            ('Benin', 'Burkina Faso', None), ('Benin', 'Niger', None),
            ('Benin', 'Nigeria', None), ('Benin', 'Togo', None),
            ('Bolivia', 'Brazil', None), ('Bolivia', 'Chile', None),
            ('Bolivia', 'Paraguay', None), ('Bolivia', 'Peru', None),
            ('Bosnia and Herzegovina', 'Croatia', None), ('Bosnia and Herzegovina', 'Montenegro', None),
            ('Bosnia and Herzegovina', 'Serbia', None), ('Brazil', 'Colombia', None),
            ('Brazil', 'Guyana', None), ('Brazil', 'Suriname', None),
            ('Brazil', 'Venezuela', None), ('Bulgaria', 'Greece', None),
            ('Bulgaria', 'Macedonia', None), ('Bulgaria', 'Romania', None),
            ('Bulgaria', 'Serbia', None), ('Burkina Faso', 'Mali', None),
            ('Burkina Faso', 'Niger', None), ('Burkina Faso', 'Togo', None),
            ('Burundi', 'Democratic Republic of the Congo', None), ('Cambodia', 'Laos', None),
            ('Cambodia', 'Thailand', None), ('Cambodia', 'Vietnam', None),
            ('Cameroon', 'Central African Republic', None), ('Cameroon', 'Chad', None),
            ('Cameroon', 'Equatorial Guinea', None), ('Cameroon', 'Nigeria', None),
            ('Cameroon', 'Republic of the Congo', None), ('Canada', 'United States', None),
            ('Central African Republic', 'Chad', None), ('Central African Republic', 'Democratic Republic of the Congo', None),
            ('Central African Republic', 'Sudan', None), ('Chad', 'Niger', None),
            ('Chad', 'Nigeria', None), ('Chad', 'Sudan', None),
            ('China', 'Bhutan', None), ('China', 'Burma', None),
            ('China', 'Hong Kong', None), ('China', 'Kazakhstan', None),
            ('China', 'Kyrgyzstan', None), ('China', 'Mongolia', None),
            ('China', 'Nepal', None), ('China', 'North Korea', None),
            ('China', 'Russia', None), ('China', 'Vietnam', None),
            ('Colombia', 'Venezuela', None), ('Costa Rica', 'Nicaragua', None),
            ("Cote d'Ivoire", 'Burkina Faso', None), ("Cote d'Ivoire", 'Guinea', None),
            ("Cote d'Ivoire", 'Mali', None), ('Cyprus', 'Akrotiri', None),
            ('Cyprus', 'Dhekelia', None), ('Czech Republic', 'Austria', None),
            ('Czech Republic', 'Germany', None), ('Czech Republic', 'Poland', None),
            ('Democratic Republic of the Congo', 'Zambia', None), ('Denmark', 'Germany', None),
            ('Djibouti', 'Eritrea', None), ('Dominican Republic', 'Haiti', None),
            ('Ecuador', 'Colombia', None), ('El Salvador', 'Honduras', None),
            ('Ethiopia', 'Djibouti', None), ('Ethiopia', 'Eritrea', None),
            ('Ethiopia', 'Kenya', None), ('Ethiopia', 'Somalia', None),
            ('Ethiopia', 'Sudan', None), ('Finland', 'Russia', None),
            ('Finland', 'Sweden', None), ('France', 'Andorra', None),
            ('France', 'Belgium', None), ('France', 'Brazil', None),
            ('France', 'Germany', None), ('France', 'Italy', None),
            ('France', 'Luxembourg', None), ('France', 'Spain', None),
            ('France', 'Suriname', None), ('France', 'Switzerland', None),
            ('Gabon', 'Cameroon', None), ('Gabon', 'Equatorial Guinea', None),
            ('Gabon', 'Republic of the Congo', None), ('Gaza Strip', 'Egypt', None),
            ('Gaza Strip', 'Israel', None), ('Ghana', 'Burkina Faso', None),
            ('Ghana', "Cote d'Ivoire", None), ('Ghana', 'Togo', None),
            ('Gibraltar', 'Spain', None), ('Guatemala', 'Belize', None),
            ('Guatemala', 'El Salvador', None), ('Guatemala', 'Honduras', None),
            ('Guatemala', 'Mexico', None), ('Guinea', 'Sierra Leone', None),
            ('Guinea-Bissau', 'Guinea', None), ('Guinea-Bissau', 'Senegal', None),
            ('Honduras', 'Nicaragua', None), ('Hungary', 'Austria', None),
            ('Hungary', 'Croatia', None), ('Hungary', 'Serbia', None),
            ('India', 'Bangladesh', None), ('India', 'Bhutan', None),
            ('India', 'Burma', None), ('India', 'China', None),
            ('India', 'Nepal', None), ('Indonesia', 'Papua New Guinea', None),
            ('Iran', 'Iraq', None), ('Ireland', 'United Kingdom', None),
            ('Israel', 'Egypt', None), ('Italy', 'Austria', None),
            ('Jordan', 'Iraq', None), ('Jordan', 'Israel', None),
            ('Jordan', 'Syria', None), ('Jordan', 'West Bank', None),
            ('Kazakhstan', 'Kyrgyzstan', None), ('Kenya', 'Somalia', None),
            ('Kenya', 'Sudan', None), ('Kenya', 'Uganda', None),
            ('Kosovo', 'Macedonia', None), ('Kosovo', 'Serbia', None),
            ('Kuwait', 'Iraq', None), ('Laos', 'Burma', None),
            ('Laos', 'China', None), ('Laos', 'Thailand', None),
            ('Laos', 'Vietnam', None), ('Latvia', 'Belarus', None),
            ('Latvia', 'Estonia', None), ('Lebanon', 'Israel', None),
            ('Lesotho', 'South Africa', None), ('Liberia', "Cote d'Ivoire", None),
            ('Liberia', 'Guinea', None), ('Liberia', 'Sierra Leone', None),
            ('Libya', 'Algeria', None), ('Libya', 'Chad', None),
            ('Libya', 'Egypt', None), ('Libya', 'Niger', None),
            ('Libya', 'Sudan', None), ('Libya', 'Tunisia', None),
            ('Liechtenstein', 'Austria', None), ('Liechtenstein', 'Switzerland', None),
            ('Lithuania', 'Belarus', None), ('Lithuania', 'Latvia', None),
            ('Lithuania', 'Poland', None), ('Lithuania', 'Russia', None),
            ('Luxembourg', 'Belgium', None), ('Luxembourg', 'Germany', None),
            ('Macau', 'China', None), ('Macedonia', 'Greece', None),
            ('Macedonia', 'Serbia', None), ('Malaysia', 'Brunei', None),
            ('Malaysia', 'Indonesia', None), ('Malaysia', 'Thailand', None),
            ('Mali', 'Algeria', None), ('Mali', 'Guinea', None),
            ('Mali', 'Niger', None), ('Mali', 'Senegal', None),
            ('Mauritania', 'Algeria', None), ('Mauritania', 'Mali', None),
            ('Mauritania', 'Senegal', None), ('Mauritania', 'Western Sahara', None),
            ('Monaco', 'France', None), ('Montenegro', 'Croatia', None),
            ('Montenegro', 'Kosovo', None), ('Montenegro', 'Serbia', None),
            ('Morocco', 'Spain', None), ('Mozambique', 'Malawi', None),
            ('Mozambique', 'Zambia', None), ('Mozambique', 'Zimbabwe', None),
            ('Namibia', 'Botswana', None), ('Namibia', 'Zambia', None),
            ('Netherlands', 'Germany', None), ('Niger', 'Algeria', None),
            ('Niger', 'Nigeria', None), ('Norway', 'Finland', None),
            ('Norway', 'Russia', None), ('Norway', 'Sweden', None),
            ('Oman', 'United Arab Emirates', None), ('Oman', 'Yemen', None),
            ('Pakistan', 'Afghanistan', None), ('Pakistan', 'China', None),
            ('Pakistan', 'India', None), ('Pakistan', 'Iran', None),
            ('Panama', 'Colombia', None), ('Panama', 'Costa Rica', None),
            ('Paraguay', 'Brazil', None), ('Peru', 'Brazil', None),
            ('Peru', 'Chile', None), ('Peru', 'Colombia', None),
            ('Peru', 'Ecuador', None), ('Poland', 'Belarus', None),
            ('Poland', 'Germany', None), ('Portugal', 'Spain', None),
            ('Republic of the Congo', 'Angola', None), ('Republic of the Congo', 'Central African Republic', None),
            ('Republic of the Congo', 'Democratic Republic of the Congo', None), ('Romania', 'Hungary', None),
            ('Romania', 'Moldova', None), ('Romania', 'Serbia', None),
            ('Russia', 'Belarus', None), ('Russia', 'Estonia', None),
            ('Russia', 'Georgia', None), ('Russia', 'Kazakhstan', None),
            ('Russia', 'Latvia', None), ('Russia', 'Mongolia', None),
            ('Russia', 'North Korea', None), ('Russia', 'Poland', None),
            ('Rwanda', 'Burundi', None), ('Rwanda', 'Democratic Republic of the Congo', None),
            ('Rwanda', 'Uganda', None), ('Saint Martin', 'Netherlands Antilles', None),
            ('San Marino', 'Italy', None), ('Saudi Arabia', 'Iraq', None),
            ('Saudi Arabia', 'Jordan', None), ('Saudi Arabia', 'Kuwait', None),
            ('Saudi Arabia', 'Oman', None), ('Saudi Arabia', 'Qatar', None),
            ('Saudi Arabia', 'United Arab Emirates', None), ('Saudi Arabia', 'Yemen', None),
            ('Senegal', 'Guinea', None), ('Serbia', 'Croatia', None),
            ('Slovakia', 'Austria', None), ('Slovakia', 'Czech Republic', None),
            ('Slovakia', 'Hungary', None), ('Slovakia', 'Poland', None),
            ('Slovakia', 'Ukraine', None), ('Slovenia', 'Austria', None),
            ('Slovenia', 'Croatia', None), ('Slovenia', 'Hungary', None),
            ('Slovenia', 'Italy', None), ('Somalia', 'Djibouti', None),
            ('South Africa', 'Botswana', None), ('South Africa', 'Mozambique', None),
            ('South Africa', 'Namibia', None), ('South Africa', 'Zimbabwe', None),
            ('South Korea', 'North Korea', None), ('Sudan', 'Democratic Republic of the Congo', None),
            ('Sudan', 'Egypt', None), ('Sudan', 'Eritrea', None),
            ('Suriname', 'Guyana', None), ('Swaziland', 'Mozambique', None),
            ('Swaziland', 'South Africa', None), ('Switzerland', 'Austria', None),
            ('Switzerland', 'Germany', None), ('Switzerland', 'Italy', None),
            ('Syria', 'Iraq', None), ('Syria', 'Israel', None),
            ('Syria', 'Lebanon', None), ('Tajikistan', 'Afghanistan', None),
            ('Tajikistan', 'China', None), ('Tajikistan', 'Kyrgyzstan', None),
            ('Tajikistan', 'Uzbekistan', None), ('Tanzania', 'Burundi', None),
            ('Tanzania', 'Democratic Republic of the Congo', None), ('Tanzania', 'Kenya', None),
            ('Tanzania', 'Malawi', None), ('Tanzania', 'Mozambique', None),
            ('Tanzania', 'Rwanda', None), ('Tanzania', 'Uganda', None),
            ('Tanzania', 'Zambia', None), ('Thailand', 'Burma', None),
            ('The Gambia', 'Senegal', None), ('Timor-Leste', 'Indonesia', None),
            ('Turkey', 'Armenia', None), ('Turkey', 'Bulgaria', None),
            ('Turkey', 'Georgia', None), ('Turkey', 'Greece', None),
            ('Turkey', 'Iran', None), ('Turkey', 'Iraq', None),
            ('Turkey', 'Syria', None), ('Turkmenistan', 'Afghanistan', None),
            ('Turkmenistan', 'Iran', None), ('Turkmenistan', 'Kazakhstan', None),
            ('Turkmenistan', 'Uzbekistan', None), ('Uganda', 'Democratic Republic of the Congo', None),
            ('Uganda', 'Sudan', None), ('Ukraine', 'Belarus', None),
            ('Ukraine', 'Hungary', None), ('Ukraine', 'Moldova', None),
            ('Ukraine', 'Poland', None), ('Ukraine', 'Romania', None),
            ('Ukraine', 'Russia', None), ('United States', 'Mexico', None),
            ('Uruguay', 'Brazil', None), ('Uzbekistan', 'Kazakhstan', None),
            ('Uzbekistan', 'Kyrgyzstan', None), ('Vatican City', 'Italy', None),
            ('Venezuela', 'Guyana', None), ('West Bank', 'Israel', None),
            ('Western Sahara', 'Algeria', None), ('Western Sahara', 'Morocco', None),
            ('Zambia', 'Malawi', None), ('Zambia', 'Zimbabwe', None),
            ('Zimbabwe', 'Botswana', None)
            ]
        gps_coordinates = {
            'Canada': [[60, 'N'], [95, 'W']],
            'Saint Martin': [[18, 'N'], [63, 'W']],
            'Sao Tome and Principe': [[1, 'N'], [7, 'E']],
            'Turkmenistan': [[40, 'N'], [60, 'E']],
            'Saint Helena': [[15, 'S'], [5, 'W']],
            'Lithuania': [[56, 'N'], [24, 'E']],
            'Cambodia': [[13, 'N'], [105, 'E']],
            'Saint Kitts and Nevis': [[17, 'N'], [62, 'W']],
            'Ethiopia': [[8, 'N'], [38, 'E']],
            'The Gambia': [[13, 'N'], [16, 'W']],
            'Aruba': [[12, 'N'], [69, 'W']],
            'Swaziland': [[26, 'S'], [31, 'E']],
            'Guinea-Bissau': [[12, 'N'], [15, 'W']],
            'Argentina': [[34, 'S'], [64, 'W']],
            'Bolivia': [[17, 'S'], [65, 'W']],
            'Bahamas, The': [[24, 'N'], [76, 'W']],
            'Spratly Islands': [[8, 'N'], [111, 'E']],
            'Ghana': [[8, 'N'], [2, 'W']],
            'Saudi Arabia': [[25, 'N'], [45, 'E']],
            'American Samoa': [[14, 'S'], [170, 'W']],
            'Cocos (Keeling) Islands': [[12, 'S'], [96, 'E']],
            'Slovenia': [[46, 'N'], [14, 'E']],
            'Guatemala': [[15, 'N'], [90, 'W']],
            'Bosnia and Herzegovina': [[44, 'N'], [18, 'E']],
            'Kuwait': [[29, 'N'], [45, 'E']],
            'Jordan': [[31, 'N'], [36, 'E']],
            'Saint Barthelemy': [[17, 'N'], [62, 'W']],
            'Ashmore and Cartier Islands': [[12, 'S'], [123, 'E']],
            'Dominica': [[15, 'N'], [61, 'W']],
            'Liberia': [[6, 'N'], [9, 'W']],
            'Maldives': [[3, 'N'], [73, 'E']],
            'Micronesia, Federated States of': [[6, 'N'], [158, 'E']],
            'Pakistan': [[30, 'N'], [70, 'E']],
            'Oman': [[21, 'N'], [57, 'E']],
            'Tanzania': [[6, 'S'], [35, 'E']],
            'Albania': [[41, 'N'], [20, 'E']],
            'Gabon': [[1, 'S'], [11, 'E']],
            'Niue': [[19, 'S'], [169, 'W']],
            'Monaco': [[43, 'N'], [7, 'E']],
            'Wallis and Futuna': [[13, 'S'], [176, 'W']],
            'New Zealand': [[41, 'S'], [174, 'E']],
            'Yemen': [[15, 'N'], [48, 'E']],
            'Jersey': [[49, 'N'], [2, 'W']],
            'Jamaica': [[18, 'N'], [77, 'W']],
            'Greenland': [[72, 'N'], [40, 'W']],
            'West Bank': [[32, 'N'], [35, 'E']],
            'Macau': [[22, 'N'], [113, 'E']],
            'Jan Mayen': [[71, 'N'], [8, 'W']],
            'United Arab Emirates': [[24, 'N'], [54, 'E']],
            'Guam': [[13, 'N'], [144, 'E']],
            'Uruguay': [[33, 'S'], [56, 'W']],
            'India': [[20, 'N'], [77, 'E']],
            'Azerbaijan': [[40, 'N'], [47, 'E']],
            'Lesotho': [[29, 'S'], [28, 'E']],
            'Saint Vincent and the Grenadines': [[13, 'N'], [61, 'W']],
            'Kenya': [[1, 'N'], [38, 'E']],
            'South Korea': [[37, 'N'], [127, 'E']],
            'Tajikistan': [[39, 'N'], [71, 'E']],
            'Turkey': [[39, 'N'], [35, 'E']],
            'Afghanistan': [[33, 'N'], [65, 'E']],
            'Paraguay': [[23, 'S'], [58, 'W']],
            'Bangladesh': [[24, 'N'], [90, 'E']],
            'Mauritania': [[20, 'N'], [12, 'W']],
            'Solomon Islands': [[8, 'S'], [159, 'E']],
            'Saint Pierre and Miquelon': [[46, 'N'], [56, 'W']],
            'Gaza Strip': [[31, 'N'], [34, 'E']],
            'San Marino': [[43, 'N'], [12, 'E']],
            'French Polynesia': [[15, 'S'], [140, 'W']],
            'France': [[46, 'N'], [2, 'E']],
            'Fiji': [[18, 'S'], [175, 'E']],
            'Rwanda': [[2, 'S'], [30, 'E']],
            'Slovakia': [[48, 'N'], [19, 'E']],
            'Somalia': [[10, 'N'], [49, 'E']],
            'Peru': [[10, 'S'], [76, 'W']],
            'Laos': [[18, 'N'], [105, 'E']],
            'Nauru': [[0, 'S'], [166, 'E']],
            'Seychelles': [[4, 'S'], [55, 'E']],
            'Norway': [[62, 'N'], [10, 'E']],
            "Cote d'Ivoire": [[8, 'N'], [5, 'W']],
            'Cook Islands': [[21, 'S'], [159, 'W']],
            'Benin': [[9, 'N'], [2, 'E']],
            'Western Sahara': [[24, 'N'], [13, 'W']],
            'Cuba': [[21, 'N'], [80, 'W']],
            'Cameroon': [[6, 'N'], [12, 'E']],
            'Montenegro': [[42, 'N'], [19, 'E']],
            'Republic of the Congo': [[1, 'S'], [15, 'E']],
            'Burkina Faso': [[13, 'N'], [2, 'W']],
            'Togo': [[8, 'N'], [1, 'E']],
            'Virgin Islands': [[18, 'N'], [64, 'W']],
            'China': [[35, 'N'], [105, 'E']],
            'Armenia': [[40, 'N'], [45, 'E']],
            'Timor-Leste': [[8, 'S'], [125, 'E']],
            'Dominican Republic': [[19, 'N'], [70, 'W']],
            'Ukraine': [[49, 'N'], [32, 'E']],
            'Bahrain': [[26, 'N'], [50, 'E']],
            'Tonga': [[20, 'S'], [175, 'W']],
            'Finland': [[64, 'N'], [26, 'E']],
            'Libya': [[25, 'N'], [17, 'E']],
            'Cayman Islands': [[19, 'N'], [80, 'W']],
            'Central African Republic': [[7, 'N'], [21, 'E']],
            'New Caledonia': [[21, 'S'], [165, 'E']],
            'Mauritius': [[20, 'S'], [57, 'E']],
            'Liechtenstein': [[47, 'N'], [9, 'E']],
            'Vietnam': [[16, 'N'], [107, 'E']],
            'British Virgin Islands': [[18, 'N'], [64, 'W']],
            'Mali': [[17, 'N'], [4, 'W']],
            'Vatican City': [[41, 'N'], [12, 'E']],
            'Russia': [[60, 'N'], [100, 'E']],
            'Bulgaria': [[43, 'N'], [25, 'E']],
            'United States': [[38, 'N'], [97, 'W']],
            'Romania': [[46, 'N'], [25, 'E']],
            'Angola': [[12, 'S'], [18, 'E']],
            'Chad': [[15, 'N'], [19, 'E']],
            'South Africa': [[29, 'S'], [24, 'E']],
            'Tokelau': [[9, 'S'], [172, 'W']],
            'Turks and Caicos Islands': [[21, 'N'], [71, 'W']],
            'South Georgia and the South Sandwich Islands': [[54, 'S'], [37, 'W']],
            'Sweden': [[62, 'N'], [15, 'E']],
            'Qatar': [[25, 'N'], [51, 'E']],
            'Malaysia': [[2, 'N'], [112, 'E']],
            'Senegal': [[14, 'N'], [14, 'W']],
            'Latvia': [[57, 'N'], [25, 'E']],
            'Clipperton Island': [[10, 'N'], [109, 'W']],
            'Uganda': [[1, 'N'], [32, 'E']],
            'Japan': [[36, 'N'], [138, 'E']],
            'Niger': [[16, 'N'], [8, 'E']],
            'Brazil': [[10, 'S'], [55, 'W']],
            'Faroe Islands': [[62, 'N'], [7, 'W']],
            'Guinea': [[11, 'N'], [10, 'W']],
            'Panama': [[9, 'N'], [80, 'W']],
            'Costa Rica': [[10, 'N'], [84, 'W']],
            'Luxembourg': [[49, 'N'], [6, 'E']],
            'Cape Verde': [[16, 'N'], [24, 'W']],
            'Andorra': [[42, 'N'], [1, 'E']],
            'Gibraltar': [[36, 'N'], [5, 'W']],
            'Ireland': [[53, 'N'], [8, 'W']],
            'Syria': [[35, 'N'], [38, 'E']],
            'Palau': [[7, 'N'], [134, 'E']],
            'Nigeria': [[10, 'N'], [8, 'E']],
            'Ecuador': [[2, 'S'], [77, 'W']],
            'Northern Mariana Islands': [[15, 'N'], [145, 'E']],
            'Brunei': [[4, 'N'], [114, 'E']],
            'Mozambique': [[18, 'S'], [35, 'E']],
            'Australia': [[27, 'S'], [133, 'E']],
            'Iran': [[32, 'N'], [53, 'E']],
            'Algeria': [[28, 'N'], [3, 'E']],
            'Svalbard': [[78, 'N'], [20, 'E']],
            'El Salvador': [[13, 'N'], [88, 'W']],
            'Tuvalu': [[8, 'S'], [178, 'E']],
            'Pitcairn Islands': [[25, 'S'], [130, 'W']],
            'Czech Republic': [[49, 'N'], [15, 'E']],
            'Marshall Islands': [[9, 'N'], [168, 'E']],
            'Chile': [[30, 'S'], [71, 'W']],
            'Puerto Rico': [[18, 'N'], [66, 'W']],
            'Belgium': [[50, 'N'], [4, 'E']],
            'Kiribati': [[1, 'N'], [173, 'E']],
            'Haiti': [[19, 'N'], [72, 'W']],
            'Belize': [[17, 'N'], [88, 'W']],
            'Hong Kong': [[22, 'N'], [114, 'E']],
            'Saint Lucia': [[13, 'N'], [60, 'W']],
            'Georgia': [[42, 'N'], [43, 'E']],
            'Mexico': [[23, 'N'], [102, 'W']],
            'Denmark': [[56, 'N'], [10, 'E']],
            'Poland': [[52, 'N'], [20, 'E']],
            'Moldova': [[47, 'N'], [29, 'E']],
            'Morocco': [[32, 'N'], [5, 'W']],
            'Namibia': [[22, 'S'], [17, 'E']],
            'Mongolia': [[46, 'N'], [105, 'E']],
            'Guernsey': [[49, 'N'], [2, 'W']],
            'Thailand': [[15, 'N'], [100, 'E']],
            'Switzerland': [[47, 'N'], [8, 'E']],
            'Grenada': [[12, 'N'], [61, 'W']],
            'Navassa Island': [[18, 'N'], [75, 'W']],
            'Isle of Man': [[54, 'N'], [4, 'W']],
            'Portugal': [[39, 'N'], [8, 'W']],
            'Estonia': [[59, 'N'], [26, 'E']],
            'Kosovo': [[42, 'N'], [21, 'E']],
            'Norfolk Island': [[29, 'S'], [167, 'E']],
            'Bouvet Island': [[54, 'S'], [3, 'E']],
            'Lebanon': [[33, 'N'], [35, 'E']],
            'Sierra Leone': [[8, 'N'], [11, 'W']],
            'Uzbekistan': [[41, 'N'], [64, 'E']],
            'Tunisia': [[34, 'N'], [9, 'E']],
            'Djibouti': [[11, 'N'], [43, 'E']],
            'Heard Island and McDonald Islands': [[53, 'S'], [72, 'E']],
            'Antigua and Barbuda': [[17, 'N'], [61, 'W']],
            'Spain': [[40, 'N'], [4, 'W']],
            'Colombia': [[4, 'N'], [72, 'W']],
            'Burundi': [[3, 'S'], [30, 'E']],
            'Taiwan': [[23, 'N'], [121, 'E']],
            'Cyprus': [[35, 'N'], [33, 'E']],
            'Barbados': [[13, 'N'], [59, 'W']],
            'Falkland Islands (Islas Malvinas)': [[51, 'S'], [59, 'W']],
            'Madagascar': [[20, 'S'], [47, 'E']],
            'Italy': [[42, 'N'], [12, 'E']],
            'Bhutan': [[27, 'N'], [90, 'E']],
            'Sudan': [[15, 'N'], [30, 'E']],
            'Vanuatu': [[16, 'S'], [167, 'E']],
            'Malta': [[35, 'N'], [14, 'E']],
            'Hungary': [[47, 'N'], [20, 'E']],
            'Democratic Republic of the Congo': [[0, 'N'], [25, 'E']],
            'Netherlands': [[52, 'N'], [5, 'E']],
            'Bermuda': [[32, 'N'], [64, 'W']],
            'Suriname': [[4, 'N'], [56, 'W']],
            'Anguilla': [[18, 'N'], [63, 'W']],
            'Venezuela': [[8, 'N'], [66, 'W']],
            'Netherlands Antilles': [[12, 'N'], [69, 'W']],
            'Israel': [[31, 'N'], [34, 'E']],
            'Paracel Islands': [[16, 'N'], [112, 'E']],
            'Wake Island': [[19, 'N'], [166, 'E']],
            'Indonesia': [[5, 'S'], [120, 'E']],
            'Iceland': [[65, 'N'], [18, 'W']],
            'Zambia': [[15, 'S'], [30, 'E']],
            'Samoa': [[13, 'S'], [172, 'W']],
            'Austria': [[47, 'N'], [13, 'E']],
            'Papua New Guinea': [[6, 'S'], [147, 'E']],
            'Malawi': [[13, 'S'], [34, 'E']],
            'Zimbabwe': [[20, 'S'], [30, 'E']],
            'Germany': [[51, 'N'], [9, 'E']],
            'Dhekelia': [[34, 'N'], [33, 'E']],
            'Kazakhstan': [[48, 'N'], [68, 'E']],
            'Philippines': [[13, 'N'], [122, 'E']],
            'Eritrea': [[15, 'N'], [39, 'E']],
            'Kyrgyzstan': [[41, 'N'], [75, 'E']],
            'Mayotte': [[12, 'S'], [45, 'E']],
            'Iraq': [[33, 'N'], [44, 'E']],
            'Montserrat': [[16, 'N'], [62, 'W']],
            'Coral Sea Islands': [[18, 'S'], [152, 'E']],
            'Macedonia': [[41, 'N'], [22, 'E']],
            'British Indian Ocean Territory': [[6, 'S'], [71, 'E']],
            'North Korea': [[40, 'N'], [127, 'E']],
            'Trinidad and Tobago': [[11, 'N'], [61, 'W']],
            'Akrotiri': [[34, 'N'], [32, 'E']],
            'Guyana': [[5, 'N'], [59, 'W']],
            'Belarus': [[53, 'N'], [28, 'E']],
            'Nepal': [[28, 'N'], [84, 'E']],
            'Burma': [[22, 'N'], [98, 'E']],
            'Honduras': [[15, 'N'], [86, 'W']],
            'Equatorial Guinea': [[2, 'N'], [10, 'E']],
            'Egypt': [[27, 'N'], [30, 'E']],
            'Nicaragua': [[13, 'N'], [85, 'W']],
            'Singapore': [[1, 'N'], [103, 'E']],
            'Serbia': [[44, 'N'], [21, 'E']],
            'Botswana': [[22, 'S'], [24, 'E']],
            'United Kingdom': [[54, 'N'], [2, 'W']],
            'Antarctica': [[90, 'S'], [0, 'E']],
            'Christmas Island': [[10, 'S'], [105, 'E']],
            'Greece': [[39, 'N'], [22, 'E']],
            'Sri Lanka': [[7, 'N'], [81, 'E']],
            'Croatia': [[45, 'N'], [15, 'E']],
            'Comoros': [[12, 'S'], [44, 'E']]
            }
        g = Graph()
        g.add_edges(edges)
        g.gps_coordinates = gps_coordinates
        g.name("World Map")
        return g



################################################################################
#   Graphs with a given degree sequence
################################################################################

    def DegreeSequence(self, deg_sequence):
        """
        Returns a graph with the given degree sequence. Raises a NetworkX
        error if the proposed degree sequence cannot be that of a graph.

        Graph returned is the one returned by the Havel-Hakimi algorithm,
        which constructs a simple graph by connecting vertices of highest
        degree to other vertices of highest degree, resorting the remaining
        vertices by degree and repeating the process. See Theorem 1.4 in
        [1].

        INPUT:


        -  ``deg_sequence`` - a list of integers with each
           entry corresponding to the degree of a different vertex.


        EXAMPLES::

            sage: G = graphs.DegreeSequence([3,3,3,3])
            sage: G.edges(labels=False)
            [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            sage: G.show()  # long time

        ::

            sage: G = graphs.DegreeSequence([3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3])
            sage: G.show()  # long time

        ::

            sage: G = graphs.DegreeSequence([4,4,4,4,4,4,4,4])
            sage: G.show()  # long time

        ::

            sage: G = graphs.DegreeSequence([1,2,3,4,3,4,3,2,3,2,1])
            sage: G.show()  # long time

        REFERENCE:

        - [1] Chartrand, G. and Lesniak, L. Graphs and Digraphs.
          Chapman and Hall/CRC, 1996.
        """
        import networkx
        return graph.Graph(networkx.havel_hakimi_graph([int(i) for i in deg_sequence]))

    def DegreeSequenceBipartite(self, s1 ,s2 ):
        r"""
        Returns a bipartite graph whose two sets have the given
        degree sequences.

        Given two different sequences of degrees `s_1` and `s_2`,
        this functions returns ( if possible ) a bipartite graph
        on sets `A` and `B` such that the vertices in `A` have
        `s_1` as their degree sequence, while `s_2` is the degree
        sequence of the vertices in `B`.

        INPUT:

        - ``s_1`` -- list of integers corresponding to the degree
          sequence of the first set.
        - ``s_2`` -- list of integers corresponding to the degree
          sequence of the second set.

        ALGORITHM:

        This function works through the computation of the matrix
        given by the Gale-Ryser theorem, which is in this case
        the adjacency matrix of the bipartite graph.

        EXAMPLES:

        If we are given as sequences ``[2,2,2,2,2]`` and ``[5,5]``
        we are given as expected the complete bipartite
        graph `K_{2,5}` ::

            sage: g = graphs.DegreeSequenceBipartite([2,2,2,2,2],[5,5])
            sage: g.is_isomorphic(graphs.CompleteBipartiteGraph(5,2))
            True

        Some sequences being incompatible if, for example, their sums
        are different, the functions raises a ``ValueError`` when no
        graph corresponding to the degree sequences exists. ::

            sage: g = graphs.DegreeSequenceBipartite([2,2,2,2,1],[5,5])
            Traceback (most recent call last):
            ...
            ValueError: There exists no bipartite graph corresponding to the given degree sequences
        """

        from sage.combinat.integer_vector import gale_ryser_theorem
        from sage.graphs.bipartite_graph import BipartiteGraph

        s1 = sorted(s1, reverse = True)
        s2 = sorted(s2, reverse = True)

        m = gale_ryser_theorem(s1,s2)

        if m is False:
            raise ValueError("There exists no bipartite graph corresponding to the given degree sequences")
        else:
            return BipartiteGraph(m)

    def DegreeSequenceConfigurationModel(self, deg_sequence, seed=None):
        """
        Returns a random pseudograph with the given degree sequence. Raises
        a NetworkX error if the proposed degree sequence cannot be that of
        a graph with multiple edges and loops.

        One requirement is that the sum of the degrees must be even, since
        every edge must be incident with two vertices.

        INPUT:


        -  ``deg_sequence`` - a list of integers with each
           entry corresponding to the expected degree of a different vertex.

        -  ``seed`` - for the random number generator.


        EXAMPLES::

            sage: G = graphs.DegreeSequenceConfigurationModel([1,1])
            sage: G.adjacency_matrix()
            [0 1]
            [1 0]

        Note: as of this writing, plotting of loops and multiple edges is
        not supported, and the output is allowed to contain both types of
        edges.

        ::

            sage: G = graphs.DegreeSequenceConfigurationModel([3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3])
            sage: G.edges(labels=False)
            [(0, 2), (0, 10), (0, 15), (1, 6), (1, 16), (1, 17), (2, 5), (2, 19), (3, 7), (3, 14), (3, 14), (4, 9), (4, 13), (4, 19), (5, 6), (5, 15), (6, 11), (7, 11), (7, 17), (8, 11), (8, 18), (8, 19), (9, 12), (9, 13), (10, 15), (10, 18), (12, 13), (12, 16), (14, 17), (16, 18)]
            sage: G.show()  # long time

        REFERENCE:

        - [1] Newman, M.E.J. The Structure and function of complex
          networks, SIAM Review vol. 45, no. 2 (2003), pp. 167-256.
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return graph.Graph(networkx.configuration_model([int(i) for i in deg_sequence], seed=seed), loops=True, multiedges=True, sparse=True)

    def DegreeSequenceTree(self, deg_sequence):
        """
        Returns a tree with the given degree sequence. Raises a NetworkX
        error if the proposed degree sequence cannot be that of a tree.

        Since every tree has one more vertex than edge, the degree sequence
        must satisfy len(deg_sequence) - sum(deg_sequence)/2 == 1.

        INPUT:


        -  ``deg_sequence`` - a list of integers with each
           entry corresponding to the expected degree of a different vertex.


        EXAMPLE::

            sage: G = graphs.DegreeSequenceTree([3,1,3,3,1,1,1,2,1])
            sage: G.show()  # long time
        """
        import networkx
        return graph.Graph(networkx.degree_sequence_tree([int(i) for i in deg_sequence]))

    def DegreeSequenceExpected(self, deg_sequence, seed=None):
        """
        Returns a random graph with expected given degree sequence. Raises
        a NetworkX error if the proposed degree sequence cannot be that of
        a graph.

        One requirement is that the sum of the degrees must be even, since
        every edge must be incident with two vertices.

        INPUT:


        -  ``deg_sequence`` - a list of integers with each
           entry corresponding to the expected degree of a different vertex.

        -  ``seed`` - for the random number generator.


        EXAMPLE::

            sage: G = graphs.DegreeSequenceExpected([1,2,3,2,3])
            sage: G.edges(labels=False)
            [(0, 2), (1, 1), (1, 3), (2, 2), (2, 4), (3, 3)]
            sage: G.show()  # long time

        REFERENCE:

        - [1] Chung, Fan and Lu, L. Connected components in random
          graphs with given expected degree
          sequences. Ann. Combinatorics (6), 2002 pp. 125-145.
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return graph.Graph(networkx.expected_degree_graph([int(i) for i in deg_sequence], seed=seed), loops=True)

###########################################################################
#   Graph Iterators
###########################################################################

    def __call__(self, vertices, property=lambda x: True, augment='edges',
        size=None, deg_seq=None, loops=False, implementation='c_graph',
        sparse=True):
        """
        Accesses the generator of isomorphism class representatives.
        Iterates over distinct, exhaustive representatives. See the docstring
        of this class for full documentation.

        EXAMPLES:

        Print graphs on 3 or less vertices::

            sage: for G in graphs(3, augment='vertices'):
            ...    print G
            Graph on 0 vertices
            Graph on 1 vertex
            Graph on 2 vertices
            Graph on 3 vertices
            Graph on 3 vertices
            Graph on 3 vertices
            Graph on 2 vertices
            Graph on 3 vertices

        For more examples, see the class level documentation, or type::

            sage: graphs? # not tested

        REFERENCE:

        - Brendan D. McKay, Isomorph-Free Exhaustive generation.
          Journal of Algorithms Volume 26, Issue 2, February 1998,
          pages 306-324.
        """
        from sage.graphs.all import Graph
        if deg_seq is not None:
            if len(deg_seq) != vertices or sum(deg_seq)%2 or sum(deg_seq) > vertices*(vertices-1):
                raise ValueError("Invalid degree sequence.")
            deg_seq = sorted(deg_seq)
            if augment == 'edges':
                property = lambda x: all([deg_seq[i] >= d for i,d in enumerate(sorted(x.degree()))])
                extra_property = lambda x: deg_seq == sorted(x.degree())
            else:
                property = lambda x: all([deg_seq[i] >= d for i,d in enumerate(sorted(x.degree() + [0]*(vertices-x.num_verts()) ))])
                extra_property = lambda x: x.num_verts() == vertices and deg_seq == sorted(x.degree())
        elif size is not None:
            extra_property = lambda x: x.size() == size
        else:
            extra_property = lambda x: True
        if augment == 'vertices':
            g = Graph(loops=loops, implementation=implementation, sparse=sparse)
            for gg in canaug_traverse_vert(g, [], vertices, property, loops=loops, implementation=implementation, sparse=sparse):
                if extra_property(gg):
                    yield gg
        elif augment == 'edges':
            g = Graph(vertices, loops=loops, implementation=implementation, sparse=sparse)
            gens = []
            for i in range(vertices-1):
                gen = range(i)
                gen.append(i+1); gen.append(i)
                gen += range(i+2, vertices)
                gens.append(gen)
            for gg in canaug_traverse_edge(g, gens, property, loops=loops, implementation=implementation, sparse=sparse):
                if extra_property(gg):
                    yield gg
        else:
            raise NotImplementedError()

    def trees(self, vertices):
        """
        Accesses the generator of trees (graphs without cycles). Iterates
        over distinct, exhaustive representatives.

        INPUT:

        -  ``vertices`` - natural number

        EXAMPLES: Sloane A000055::

            sage: for i in range(0, 7):
            ...    print len(list(graphs.trees(i)))
            1
            1
            1
            1
            2
            3
            6
            sage: for i in range(7, 10):            # long time
            ...    print len(list(graphs.trees(i))) # long time
            11
            23
            47
        """
        from trees import TreeIterator
        return iter(TreeIterator(vertices))

    def nauty_geng(self, options=""):
        r"""
        Calls the geng program in the optional nauty spkg to generate
        graphs. The options argument is passed straight to nauty.

        INPUT:


        -  ``options`` - a string passed to the command line of
           geng. You *must* pass the number of vertices you desire.


        EXAMPLES::

            sage: graph_list = graphs.nauty_geng("-q 3") # requires the optional nauty package
            sage: len(graph_list) # requires the optional nauty package
            4
        """
        import os
        from sage.misc.package import is_package_installed
        if not is_package_installed("nauty"):
            raise TypeError, "the optional nauty package is not installed"
        return [graph.Graph(g) for g in os.popen("nauty-geng %s"%(options) ).read().split()]

    def cospectral_graphs(self, vertices, matrix_function=lambda g: g.adjacency_matrix(), graphs=None):
        """
        Find all sets of graphs on ``vertices`` vertices (with
        possible restrictions) which are cospectral with respect to a
        constructed matrix.

        INPUT:

        - ``vertices`` - The number of vertices in the graphs to be tested

        - ``matrix_function`` - A function taking a graph and giving back
          a matrix.  This defaults to the adjacency matrix.  The spectra
          examined are the spectra of these matrices.

        - ``graphs`` - One of three things:

           - ``None`` (default) - test all graphs having ``vertices``
             vertices

           - a function taking a graph and returning ``True`` or ``False``
             - test only the graphs on ``vertices`` vertices for which
             the function returns ``True``

           - a list of graphs (or other iterable object) - these graphs
             are tested for cospectral sets.  In this case,
             ``vertices`` is ignored.

        OUTPUT:

           A list of lists of graphs.  Each sublist will be a list of
           cospectral graphs.


        EXAMPLES::

            sage: g=graphs.cospectral_graphs(5)
            sage: sorted(sorted(g.graph6_string() for g in glist) for glist in g)
            [['Dr?', 'Ds_']]
            sage: g[0][1].am().charpoly()==g[0][1].am().charpoly()
            True

        There are two sets of cospectral graphs on six vertices with no isolated vertices::

            sage: g=graphs.cospectral_graphs(6, graphs=lambda x: min(x.degree())>0)
            sage: sorted(sorted(g.graph6_string() for g in glist) for glist in g)
            [['Ep__', 'Er?G'], ['ExGg', 'ExoG']]
            sage: g[0][1].am().charpoly()==g[0][1].am().charpoly()
            True
            sage: g[1][1].am().charpoly()==g[1][1].am().charpoly()
            True

        There is one pair of cospectral trees on eight vertices::

            sage: g=graphs.cospectral_graphs(6, graphs=graphs.trees(8))
            sage: sorted(sorted(g.graph6_string() for g in glist) for glist in g)
            [['GiPC?C', 'GiQCC?']]
            sage: g[0][1].am().charpoly()==g[0][1].am().charpoly()
            True

        There are two sets of cospectral graphs (with respect to the
        Laplacian matrix) on six vertices::

            sage: g=graphs.cospectral_graphs(6, matrix_function=lambda g: g.laplacian_matrix())
            sage: sorted(sorted(g.graph6_string() for g in glist) for glist in g)
            [['Edq_', 'ErcG'], ['Exoo', 'EzcG']]
            sage: g[0][1].laplacian_matrix().charpoly()==g[0][1].laplacian_matrix().charpoly()
            True
            sage: g[1][1].laplacian_matrix().charpoly()==g[1][1].laplacian_matrix().charpoly()
            True

        To find cospectral graphs with respect to the normalized
        Laplacian, assuming the graphs do not have an isolated vertex, it
        is enough to check the spectrum of the matrix `D^{-1}A`, where D
        is the diagonal matrix of vertex degrees, and A is the adjacency
        matrix.  We find two such cospectral graphs (for the normalized
        Laplacian) on five vertices::

            sage: def DinverseA(g):
            ...     A=g.adjacency_matrix().change_ring(QQ)
            ...     for i in range(g.order()):
            ...         A.rescale_row(i, 1/len(A.nonzero_positions_in_row(i)))
            ...     return A
            sage: g=graphs.cospectral_graphs(5, matrix_function=DinverseA, graphs=lambda g: min(g.degree())>0)
            sage: sorted(sorted(g.graph6_string() for g in glist) for glist in g)
            [['Dlg', 'Ds_']]
            sage: g[0][1].laplacian_matrix(normalized=True).charpoly()==g[0][1].laplacian_matrix(normalized=True).charpoly()
            True
        """
        from sage.graphs.all import graphs as graph_gen
        if graphs is None:
            graph_list=graph_gen(vertices)
        elif callable(graphs):
            graph_list=iter(g for g in graph_gen(vertices) if graphs(g))
        else:
            graph_list=iter(graphs)

        from collections import defaultdict
        charpolys=defaultdict(list)
        for g in graph_list:
            cp=matrix_function(g).charpoly()
            charpolys[cp].append(g)

        cospectral_graphs=[]
        for cp,g_list in charpolys.items():
            if len(g_list)>1:
                cospectral_graphs.append(g_list)

        return cospectral_graphs


def canaug_traverse_vert(g, aut_gens, max_verts, property, dig=False, loops=False, implementation='c_graph', sparse=True):
    """
    Main function for exhaustive generation. Recursive traversal of a
    canonically generated tree of isomorph free (di)graphs satisfying a
    given property.

    INPUT:


    -  ``g`` - current position on the tree.

    -  ``aut_gens`` - list of generators of Aut(g), in
       list notation.

    -  ``max_verts`` - when to retreat.

    -  ``property`` - check before traversing below g.

    -  ``deg_seq`` - specify a degree sequence to try to
       obtain.


    EXAMPLES::

        sage: from sage.graphs.graph_generators import canaug_traverse_vert
        sage: list(canaug_traverse_vert(Graph(), [], 3, lambda x: True))
        [Graph on 0 vertices, ... Graph on 3 vertices]

    The best way to access this function is through the graphs()
    iterator:

    Print graphs on 3 or less vertices.

    ::

        sage: for G in graphs(3, augment='vertices'):
        ...    print G
        ...
        Graph on 0 vertices
        Graph on 1 vertex
        Graph on 2 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 2 vertices
        Graph on 3 vertices

    Print digraphs on 2 or less vertices.

    ::

        sage: for D in digraphs(2, augment='vertices'):
        ...    print D
        ...
        Digraph on 0 vertices
        Digraph on 1 vertex
        Digraph on 2 vertices
        Digraph on 2 vertices
        Digraph on 2 vertices
    """
    from sage.graphs.generic_graph_pyx import binary
    from sage.groups.perm_gps.partn_ref.refinement_graphs import search_tree


    if not property(g):
        return
    yield g

    n = g.order()
    if n < max_verts:

        # build a list representing C(g) - the vertex to be added
        # is at the end, so only specify which edges...
        # in the case of graphs, there are n possibilities,
        # and in the case of digraphs, there are 2*n.
        if dig:
            possibilities = 2*n
        else:
            possibilities = n
        num_roots = 2**possibilities
        children = [-1]*num_roots

        # union-find C(g) under Aut(g)
        for gen in aut_gens:
            for i in xrange(len(children)):
                k = 0
                for j in xrange(possibilities):
                    if (1 << j)&i:
                        if dig and j >= n:
                            k += (1 << (gen[j-n]+n))
                        else:
                            k += (1 << gen[j])
                while children[k] != -1:
                    k = children[k]
                while children[i] != -1:
                    i = children[i]
                if i != k:
                    # union i & k
                    smaller, larger = sorted([i,k])
                    children[larger] = smaller
                    num_roots -= 1

        # find representatives of orbits of C(g)
        roots = []
        found_roots = 0
        i = 0
        while found_roots < num_roots:
            if children[i] == -1:
                found_roots += 1
                roots.append(i)
            i += 1
        for i in roots:
            # construct a z for each number in roots...
            z = g.copy(implementation=implementation, sparse=sparse)
            z.add_vertex(n)
            edges = []
            if dig:
                index = 0
                while index < possibilities/2:
                    if (1 << index)&i:
                        edges.append((index,n))
                    index += 1
                while index < possibilities:
                    if (1 << index)&i:
                        edges.append((n,index-n))
                    index += 1
            else:
                index = 0
                while (1 << index) <= i:
                    if (1 << index)&i:
                        edges.append((index,n))
                    index += 1
            z.add_edges(edges)
            z_s = []
            if property(z):
                z_s.append(z)
            if loops:
                z = z.copy(implementation=implementation, sparse=sparse)
                z.add_edge((n,n))
                if property(z):
                    z_s.append(z)
            for z in z_s:
                z_aut_gens, _, canonical_relabeling = search_tree(z, [z.vertices()], certify=True, dig=(dig or loops))
                cut_vert = 0
                while canonical_relabeling[cut_vert] != n:
                    cut_vert += 1
                sub_verts = [v for v in z if v != cut_vert]
                m_z = z.subgraph(sub_verts)

                if m_z == g:
                    for a in canaug_traverse_vert(z, z_aut_gens, max_verts, property, dig=dig, loops=loops, implementation=implementation, sparse=sparse):
                        yield a
                else:
                    for possibility in check_aut(z_aut_gens, cut_vert, n):
                        if m_z.relabel(possibility, inplace=False) == g:
                            for a in canaug_traverse_vert(z, z_aut_gens, max_verts, property, dig=dig, loops=loops, implementation=implementation, sparse=sparse):
                                yield a
                            break

def check_aut(aut_gens, cut_vert, n):
    """
    Helper function for exhaustive generation.

    At the start, check_aut is given a set of generators for the
    automorphism group, aut_gens. We already know we are looking for
    an element of the auto- morphism group that sends cut_vert to n,
    and check_aut generates these for the canaug_traverse function.

    EXAMPLE: Note that the last two entries indicate that none of the
    automorphism group has yet been searched - we are starting at the
    identity [0, 1, 2, 3] and so far that is all we have seen. We
    return automorphisms mapping 2 to 3.

    ::

        sage: from sage.graphs.graph_generators import check_aut
        sage: list( check_aut( [ [0, 3, 2, 1], [1, 0, 3, 2], [2, 1, 0, 3] ], 2, 3))
        [[1, 0, 3, 2], [1, 2, 3, 0]]
    """
    from copy import copy
    perm = range(n+1)
    seen_perms = [perm]
    unchecked_perms = [perm]
    while len(unchecked_perms) != 0:
        perm = unchecked_perms.pop(0)
        for gen in aut_gens:
            new_perm = copy(perm)
            for i in xrange(len(perm)):
                new_perm[i] = gen[perm[i]]
            if new_perm not in seen_perms:
                seen_perms.append(new_perm)
                unchecked_perms.append(new_perm)
                if new_perm[cut_vert] == n:
                    yield new_perm

def canaug_traverse_edge(g, aut_gens, property, dig=False, loops=False, implementation='c_graph', sparse=True):
    """
    Main function for exhaustive generation. Recursive traversal of a
    canonically generated tree of isomorph free graphs satisfying a
    given property.

    INPUT:


    -  ``g`` - current position on the tree.

    -  ``aut_gens`` - list of generators of Aut(g), in
       list notation.

    -  ``property`` - check before traversing below g.


    EXAMPLES::

        sage: from sage.graphs.graph_generators import canaug_traverse_edge
        sage: G = Graph(3)
        sage: list(canaug_traverse_edge(G, [], lambda x: True))
        [Graph on 3 vertices, ... Graph on 3 vertices]

    The best way to access this function is through the graphs()
    iterator:

    Print graphs on 3 or less vertices.

    ::

        sage: for G in graphs(3):
        ...    print G
        ...
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices

    Print digraphs on 3 or less vertices.

    ::

        sage: for G in digraphs(3):
        ...    print G
        ...
        Digraph on 3 vertices
        Digraph on 3 vertices
        ...
        Digraph on 3 vertices
        Digraph on 3 vertices
    """
    from sage.graphs.generic_graph_pyx import binary
    from sage.groups.perm_gps.partn_ref.refinement_graphs import search_tree
    if not property(g):
        return
    yield g
    n = g.order()
    if dig:
        max_size = n*(n-1)
    else:
        max_size = (n*(n-1))>>1 # >> 1 is just / 2 (this is n choose 2)
    if loops: max_size += n
    if g.size() < max_size:
        # build a list representing C(g) - the edge to be added
        # is one of max_size choices
        if dig:
            children = [[(j,i) for i in xrange(n)] for j in xrange(n)]
        else:
            children = [[(j,i) for i in xrange(j)] for j in xrange(n)]
        # union-find C(g) under Aut(g)
        orbits = range(n)
        for gen in aut_gens:
            for iii in xrange(n):
                if orbits[gen[iii]] != orbits[iii]:
                    temp = orbits[gen[iii]]
                    for jjj in xrange(n):
                        if orbits[jjj] == temp:
                            orbits[jjj] = orbits[iii]
                if dig:
                    jjj_range = range(iii) + range(iii+1, n)
                else:
                    jjj_range = xrange(iii) # iii > jjj
                for jjj in jjj_range:
                    i, j = iii, jjj
                    if dig:
                        x, y = gen[i], gen[j]
                    else:
                        y, x = sorted([gen[i], gen[j]])
                    if children[i][j] != children[x][y]:
                        x_val, y_val = x, y
                        i_val, j_val = i, j
                        if dig:
                            while (x_val, y_val) != children[x_val][y_val]:
                                x_val, y_val = children[x_val][y_val]
                            while (i_val, j_val) != children[i_val][j_val]:
                                i_val, j_val = children[i_val][j_val]
                        else:
                            while (x_val, y_val) != children[x_val][y_val]:
                                y_val, x_val = sorted(children[x_val][y_val])
                            while (i_val, j_val) != children[i_val][j_val]:
                                j_val, i_val = sorted(children[i_val][j_val])
                        while (x, y) != (x_val, y_val):
                            xx, yy = x, y
                            x, y = children[x][y]
                            children[xx][yy] = (x_val, y_val)
                        while (i, j) != (i_val, j_val):
                            ii, jj = i, j
                            i, j = children[i][j]
                            children[ii][jj] = (i_val, j_val)
                        if x < i:
                            children[i][j] = (x, y)
                        elif x > i:
                            children[x][y] = (i, j)
                        elif y < j:
                            children[i][j] = (x, y)
                        elif y > j:
                            children[x][y] = (i, j)
                        else:
                            continue
        # find representatives of orbits of C(g)
        roots = []
        for i in range(n):
            if dig:
                j_range = range(i) + range(i+1, n)
            else:
                j_range = range(i)
            for j in j_range:
                if children[i][j] == (i, j):
                    roots.append((i,j))
        if loops:
            seen = []
            for i in xrange(n):
                if orbits[i] not in seen:
                    roots.append((i,i))
                    seen.append(orbits[i])
        for i, j in roots:
            if g.has_edge(i, j):
                continue
            # construct a z for each edge in roots...
            z = g.copy(implementation=implementation, sparse=sparse)
            z.add_edge(i, j)
            if not property(z):
                continue
            z_aut_gens, _, canonical_relabeling = search_tree(z, [z.vertices()], certify=True, dig=(dig or loops))
            relabel_inverse = [0]*n
            for ii in xrange(n):
                relabel_inverse[canonical_relabeling[ii]] = ii
            z_can = z.relabel(canonical_relabeling, inplace=False)
            cut_edge_can = z_can.edges(labels=False, sort=True)[-1]
            cut_edge = [relabel_inverse[cut_edge_can[0]], relabel_inverse[cut_edge_can[1]]]
            if dig:
                cut_edge = tuple(cut_edge)
            else:
                cut_edge = tuple(sorted(cut_edge))

            from copy import copy
            m_z = copy(z)
            m_z.delete_edge(cut_edge)
            if m_z == g:
                for a in canaug_traverse_edge(z, z_aut_gens, property, dig=dig, loops=loops, implementation=implementation, sparse=sparse):
                    yield a
            else:
                for possibility in check_aut_edge(z_aut_gens, cut_edge, i, j, n, dig=dig):
                    if m_z.relabel(possibility, inplace=False) == g:
                        for a in canaug_traverse_edge(z, z_aut_gens, property, dig=dig, loops=loops, implementation=implementation, sparse=sparse):
                            yield a
                        break

def check_aut_edge(aut_gens, cut_edge, i, j, n, dig=False):
    """
    Helper function for exhaustive generation.

    At the start, check_aut_edge is given a set of generators for the
    automorphism group, aut_gens. We already know we are looking for
    an element of the auto- morphism group that sends cut_edge to {i,
    j}, and check_aut generates these for the canaug_traverse
    function.

    EXAMPLE: Note that the last two entries indicate that none of the
    automorphism group has yet been searched - we are starting at the
    identity [0, 1, 2, 3] and so far that is all we have seen. We
    return automorphisms mapping 2 to 3.

    ::

        sage: from sage.graphs.graph_generators import check_aut
        sage: list( check_aut( [ [0, 3, 2, 1], [1, 0, 3, 2], [2, 1, 0, 3] ], 2, 3))
        [[1, 0, 3, 2], [1, 2, 3, 0]]
    """
    from copy import copy
    perm = range(n)
    seen_perms = [perm]
    unchecked_perms = [perm]
    while len(unchecked_perms) != 0:
        perm = unchecked_perms.pop(0)
        for gen in aut_gens:
            new_perm = copy(perm)
            for ii in xrange(n):
                new_perm[ii] = gen[perm[ii]]
            if new_perm not in seen_perms:
                seen_perms.append(new_perm)
                unchecked_perms.append(new_perm)
                if new_perm[cut_edge[0]] == i and new_perm[cut_edge[1]] == j:
                    yield new_perm
                if not dig and new_perm[cut_edge[0]] == j and new_perm[cut_edge[1]] == i:
                    yield new_perm


# Easy access to the graph generators from the command line:
graphs = GraphGenerators()




