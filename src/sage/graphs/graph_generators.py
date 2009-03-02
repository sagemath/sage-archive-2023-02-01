r"""
A collection of constructors of common graphs.

USE:

To see a list of all graph constructors, type "graphs." and then
press the tab key. The documentation for each constructor includes
information about each graph, which provides a useful reference.

PLOTTING:

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
below, particularly for CycleGraph, StarGraph, WheelGraph,
CompleteGraph and CompleteBipartiteGraph.

ORGANIZATION:

The constructors available in this database are
organized as follows::

    Basic Structures:
        - BarbellGraph
        - BullGraph
        - CircularLadderGraph
        - ClawGraph
        - CycleGraph
        - DiamondGraph
        - EmptyGraph
        - Grid2dGraph
        - GridGraph
        - HouseGraph
        - HouseXGraph
        - KrackhardtKiteGraph
        - LadderGraph
        - LollipopGraph
        - PathGraph
        - StarGraph
        - WheelGraph
    Platonic Solids:
        - TetrahedralGraph
        - HexahedralGraph
        - OctahedralGraph
        - IcosahedralGraph
        - DodecahedralGraph
    Named Graphs:
        - ChvatalGraph
        - DesarguesGraph
        - FlowerSnark
        - FruchtGraph
        - HeawoodGraph
        - HoffmanSingletonGraph
        - MoebiusKantorGraph
        - Pappus Graph
        - PetersenGraph
        - ThomsenGraph
    Families of Graphs:
        - CirculantGraph
        - CompleteGraph
        - CompleteBipartiteGraph
        - CubeGraph
        - BalancedTree
        - LCFGraph
    Pseudofractal Graphs:
        - DorogovtsevGoltsevMendesGraph
    Random Graphs:
        - RandomGNP
        - RandomBarabasiAlbert
        - RandomGNM
        - RandomNewmanWattsStrogatz
        - RandomHolmeKim
        - RandomLobster
        - RandomTreePowerlaw
        - RandomRegular
        - RandomShell
    Random Directed Graphs:
        - RandomDirectedGN
        - RandomDirectedGNC
        - RandomDirectedGNR
    Graphs with a given degree sequence:
        - DegreeSequence
        - DegreeSequenceConfigurationModel
        - DegreeSequenceTree
        - DegreeSequenceExpected


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
"""

################################################################################
#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
################################################################################

import graph
from   math import sin, cos, pi
from sage.misc.randstate import current_randstate

class GraphGenerators():
    r"""
    A class consisting of constructors for several common graphs, as
    well as orderly generation of isomorphism class representatives.

    A list of all graphs and graph structures (other than iso. class
    rep's) in this database is available via tab completion. Type
    "graphs." and then hit tab to see which graphs are available.

    The docstrings include educational information about each named
    graph with the hopes that this class can be used as a reference.

    For all the constructors in this class (except the octahedral,
    dodecahedral, random and empty graphs), the position dictionary is
    filled to override the spring-layout algorithm.

    The constructors currently in this class include::

                Basic Structures:
                    - BarbellGraph
                    - BullGraph
                    - CircularLadderGraph
                    - ClawGraph
                    - CycleGraph
                    - DiamondGraph
                    - EmptyGraph
                    - Grid2dGraph
                    - GridGraph
                    - HouseGraph
                    - HouseXGraph
                    - KrackhardtKiteGraph
                    - LadderGraph
                    - LollipopGraph
                    - PathGraph
                    - StarGraph
                    - WheelGraph
                Platonic Solids:
                    - TetrahedralGraph
                    - HexahedralGraph
                    - OctahedralGraph
                    - IcosahedralGraph
                    - DodecahedralGraph
                Named Graphs:
                    - ChvatalGraph
                    - DesarguesGraph
                    - FlowerSnark
                    - FruchtGraph
                    - HeawoodGraph
                    - HoffmanSingletonGraph
                    - MoebiusKantorGraph
                    - PappusGraph
                    - PetersenGraph
                    - ThomsenGraph
                Families of Graphs:
                    - CirculantGraph
                    - CompleteGraph
                    - CompleteBipartiteGraph
                    - CubeGraph
                    - BalancedTree
                    - LCFGraph
                Pseudofractal Graphs:
                    - DorogovtsevGoltsevMendesGraph
                Random Graphs:
                    - RandomGNP
                    - RandomBarabasiAlbert
                    - RandomGNM
                    - RandomNewmanWattsStrogatz
                    - RandomHolmeKim
                    - RandomLobster
                    - RandomTreePowerlaw
                    - RandomRegular
                    - RandomShell
                Graphs with a given degree sequence:
                    - DegreeSequence
                    - DegreeSequenceConfigurationModel
                    - DegreeSequenceTree
                    - DegreeSequenceExpected


    ORDERLY GENERATION: graphs(vertices, property=lambda x: True,
    augment='edges', size=None)

    This syntax accesses the generator of isomorphism class
    representatives. Iterates over distinct, exhaustive
    representatives.

    INPUT:


    -  ``vertices`` - natural number

    -  ``property`` - any property to be tested on graphs
       before generation. (Ignored if deg_seq is specified.)

    -  ``augment`` - choices:

    -  ``'vertices'`` - augments by adding a vertex, and
       edges incident to that vertex. In this case, all graphs on up to
       n=vertices are generated. If for any graph G satisfying the
       property, every subgraph, obtained from G by deleting one vertex
       and only edges incident to that vertex, satisfies the property,
       then this will generate all graphs with that property. If this does
       not hold, then all the graphs generated will satisfy the property,
       but there will be some missing.

    -  ``'edges'`` - augments a fixed number of vertices by
       adding one edge In this case, all graphs on exactly n=vertices are
       generated. If for any graph G satisfying the property, every
       subgraph, obtained from G by deleting one edge but not the vertices
       incident to that edge, satisfies the property, then this will
       generate all graphs with that property. If this does not hold, then
       all the graphs generated will satisfy the property, but there will
       be some missing.

    -  ``deg_seq`` - a sequence of non-negative integers,
       or None. If specified, the generated graphs will have these
       integers for degrees. In this case property and size are both
       ignored.

    -  ``loops`` - whether to allow loops in the graph or
       not.

    -  ``implementation`` - which underlying implementation to use (see Graph?)

    -  ``sparse`` - ignored if implementation is not ``c_graph``

    EXAMPLES: Print graphs on 3 or less vertices.

    ::

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

        sage: L = list(graphs(6,augment='vertices',loops=True))               # long time
        sage: for i in [0..6]: print i, len([g for g in L if g.order() == i]) # long time
        0 1
        1 2
        2 6
        3 20
        4 90
        5 544
        6 5096

    Generate all graphs with a specified degree sequence: (see
    http://www.research.att.com/~njas/sequences/A002851)

    ::

        sage: for i in [4,6,8]:
        ...    print i, len([g for g in graphs(i,deg_seq=[3]*i) if g.is_connected()])
        4 1
        6 2
        8 5
        sage: for i in [4,6,8]:                                                                          # long time
        ...    print i, len([g for g in graphs(i,augment='vertices',deg_seq=[3]*i) if g.is_connected()]) # long time
        4 1
        6 2
        8 5

    ::

        sage: print 10, len([g for g in graphs(10,deg_seq=[3]*10) if g.is_connected()]) # not tested
        10 19

    REFERENCE:

    - Brendan D. McKay, Isomorph-Free Exhaustive generation.  Journal
      of Algorithms Volume 26, Issue 2, February 1998, pages 306-324.
    """

################################################################################
#   Basic Structures
################################################################################

    def BarbellGraph(self, n1, n2):
        """
        Returns a barbell graph with 2\*n1 + n2 nodes. n1 must be greater
        than or equal to 2.

        A barbell graph is a basic structure that consists of a path graph
        of order n2 connecting two complete graphs of order n1 each.

        This constructor depends on NetworkX numeric labels. In this case,
        the (n1)th node connects to the path graph from one complete graph
        and the (n1+n2+1)th node connects to the path graph from the other
        complete graph.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, each barbell
        graph will be displayed with the two complete graphs in the
        lower-left and upper-right corners, with the path graph connecting
        diagonally between the two. Thus the (n1)th node will be drawn at a
        45 degree angle from the horizontal right center of the first
        complete graph, and the (n1+n2+1)th node will be drawn 45 degrees
        below the left horizontal center of the second complete graph.

        EXAMPLES: Construct and show a barbell graph Bar = 4, Bells = 9

        ::

            sage: g = graphs.BarbellGraph(9,4)
            sage: g.show() # long time

        Create several barbell graphs in a Sage graphics array

        ::

            sage: g = []
            sage: j = []
            sage: for i in range(6):
            ...    k = graphs.BarbellGraph(i+2,4)
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
            pos_dict[j] = [x,y]
        for i in range(n1+n2)[n1:]:
            x = float(i - n1 - n2/2 + 1)
            y = float(i - n1 - n2/2 + 1)
            pos_dict[i] = [x,y]
        for i in range(2*n1+n2)[n1+n2:]:
            x = float(cos((5*pi/4) + ((2*pi)/n1)*(i-n1-n2)) + n2/2 + 2)
            y = float(sin((5*pi/4) + ((2*pi)/n1)*(i-n1-n2)) + n2/2 + 2)
            pos_dict[i] = [x,y]

        import networkx
        G = networkx.barbell_graph(n1,n2)
        return graph.Graph(G, pos=pos_dict, name="Barbell graph")

    def BullGraph(self):
        """
        Returns a bull graph with 5 nodes.

        A bull graph is named for its shape. It's a triangle with horns.

        This constructor depends on NetworkX numeric labeling.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the bull graph
        is drawn as a triangle with the first node (0) on the bottom. The
        second and third nodes (1 and 2) complete the triangle. Node 3 is
        the horn connected to 1 and node 4 is the horn connected to node
        2.

        EXAMPLES: Construct and show a bull graph

        ::

            sage: g = graphs.BullGraph()
            sage: g.show() # long time
        """
        pos_dict = {0:[0,0],1:[-1,1],2:[1,1],3:[-2,2],4:[2,2]}
        import networkx
        G = networkx.bull_graph()
        return graph.Graph(G, pos=pos_dict, name="Bull Graph")


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
        for i in range(2*n)[n:]:
            x = float(2*(cos((pi/2) + ((2*pi)/n)*(i-n))))
            y = float(2*(sin((pi/2) + ((2*pi)/n)*(i-n))))
            pos_dict[i] = [x,y]
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
        pos_dict = {0:[0,1],1:[-1,0],2:[0,0],3:[1,0]}
        import networkx
        G = networkx.complete_bipartite_graph(1,3)
        return graph.Graph(G, pos=pos_dict, name="Claw graph")

    def CycleGraph(self, n):
        r"""
        Returns a cycle graph with n nodes.

        A cycle graph is a basic structure which is also typically called
        an n-gon.

        This constructor is dependant on vertices numbered 0 through n-1 in
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
            pos_dict[i] = [x,y]
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
        pos_dict = {0:[0,1],1:[-1,0],2:[1,0],3:[0,-1]}
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

    def Grid2dGraph(self, n1, n2):
        """
        Returns a 2-dimensional grid graph with n1\*n2 nodes (n1 rows and
        n2 columns).

        A 2d grid graph resembles a 2 dimensional grid. All inner nodes are
        connected to their 4 neighbors. Outer (non-corner) nodes are
        connected to their 3 neighbors. Corner nodes are connected to their
        2 neighbors.

        This constructor depends on NetworkX numeric labels.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, nodes are
        labelled in (row, column) pairs with (0, 0) in the top left corner.
        Edges will always be horizontal and vertical - another advantage of
        filling the position dictionary.

        EXAMPLES: Construct and show a grid 2d graph Rows = 5, Columns = 7

        ::

            sage: g = graphs.Grid2dGraph(5,7)
            sage: g.show() # long time
        """
        pos_dict = {}
        for i in range(n1):
            y = -i
            for j in range(n2):
                x = j
                pos_dict[i,j] = [x,y]
        import networkx
        G = networkx.grid_2d_graph(n1,n2)
        return graph.Graph(G, pos=pos_dict, name="2D Grid Graph", implementation='networkx')

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
        return graph.Graph(G, name="Grid Graph for %s"%dim, implementation='networkx')

    def HouseGraph(self):
        """
        Returns a house graph with 5 nodes.

        A house graph is named for its shape. It is a triange (roof) over a
        square (walls).

        This constructor depends on NetworkX numeric labeling.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the house
        graph is drawn with the first node in the lower-left corner of the
        house, the second in the lower-right corner of the house. The third
        node is in the upper-left corner connecting the roof to the wall,
        and the fourth is in the upper-right corner connecting the roof to
        the walll. The fifth node is the top of the roof, connected only to
        the third and fourth.

        EXAMPLES: Construct and show a house graph

        ::

            sage: g = graphs.HouseGraph()
            sage: g.show() # long time
        """
        pos_dict = {0:[-1,0],1:[1,0],2:[-1,1],3:[1,1],4:[0,2]}
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
        the walll. The fifth node is the top of the roof, connected only to
        the third and fourth.

        EXAMPLES: Construct and show a house X graph

        ::

            sage: g = graphs.HouseXGraph()
            sage: g.show() # long time
        """
        pos_dict = {0:[-1,0],1:[1,0],2:[-1,1],3:[1,1],4:[0,2]}
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
        clique (i.e.: Degree Centrality). The eigth (7) node is where the
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
        pos_dict = {0:[-1,4],1:[1,4],2:[-2,3],3:[0,3],4:[2,3],5:[-1,2],6:[1,2],7:[0,1],8:[0,0],9:[0,-1]}
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
            pos_dict[i] = [i,1]
        for i in range(2*n)[n:]:
            x = i - n
            pos_dict[i] = [x,0]
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
            pos_dict[j] = [x,y]
        for i in range(n1+n2)[n1:]:
            x = float(i - n1 - n2/2 + 1)
            y = float(i - n1 - n2/2 + 1)
            pos_dict[i] = [x,y]

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
                pos_dict[i] = [x,y]
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
                    pos_dict[counter] = [x,y]
                    counter += 1
                if lr: lr = False
                else: lr = True
            y = -rows
            for j in range(rem): # last row
                if lr:
                    x = j
                else:
                    x = 9 - j
                pos_dict[counter] = [x,y]
                counter += 1

        import networkx
        G = networkx.path_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Path Graph")

    def StarGraph(self, n):
        """
        Returns a star graph with n+1 nodes.

        A Star graph is a basic structure where one node is connected to
        all other nodes.

        This constructor is dependant on NetworkX numeric labels.

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
        pos_dict[0] = [0,0]
        for i in range(n+1)[1:]:
            x = float(cos((pi/2) + ((2*pi)/n)*(i-1)))
            y = float(sin((pi/2) + ((2*pi)/n)*(i-1)))
            pos_dict[i] = [x,y]
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
        pos_dict[0] = [0,0]
        for i in range(n)[1:]:
            x = float(cos((pi/2) + ((2*pi)/(n-1))*(i-1)))
            y = float(sin((pi/2) + ((2*pi)/(n-1))*(i-1)))
            pos_dict[i] = [x,y]
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

        EXAMPLES: Construct and show a Dodecahdedral graph

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

################################################################################
#   Named Graphs
################################################################################

    def ChvatalGraph(self):
        """
        Returns the Chvatal graph.

        The Chvatal graph has 12 vertices. It is a 4-regular, 4-chromatic
        graph. It is one of the few known graphs to satisfy Grunbaum's
        conjecture that for every m 1, n 2, there is an m-regular,
        m-chromatic graph of girth at least n.

        EXAMPLE::

            sage: G = graphs.ChvatalGraph()
            sage: G.degree()
            [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
        """
        import networkx
        pos_dict = {}
        for i in range(10)[5:]:
            x = float(cos((pi/2) + ((2*pi)/5)*i))
            y = float(sin((pi/2) + ((2*pi)/5)*i))
            pos_dict[i] = [x,y]
        for i in range(5):
            x = float(2*(cos((pi/2) + ((2*pi)/5)*(i-5))))
            y = float(2*(sin((pi/2) + ((2*pi)/5)*(i-5))))
            pos_dict[i] = [x,y]
        pos_dict[10] = [.5,0]
        pos_dict[11] = [-.5,0]

        return graph.Graph(networkx.chvatal_graph(), pos=pos_dict, name="Chvatal Graph")

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
        pos_dict = {}
        for i in range(10):
            x = float(cos(pi/2 + ((2*pi)/10)*i))
            y = float(sin(pi/2 + ((2*pi)/10)*i))
            pos_dict[i] = [x,y]
        for i in range(20)[10:]:
            x = float(0.5*cos(pi/2 + ((2*pi)/10)*i))
            y = float(0.5*sin(pi/2 + ((2*pi)/10)*i))
            pos_dict[i] = [x,y]
        G = graph.Graph({0:[1,9,10], 1:[2,11], 2:[3,12], 3:[4,13], 4:[5,14],\
                   5:[6,15], 6:[7,16], 7:[8,17], 8:[9,18], 9:[19], 10:[13,17],\
                   11:[14,18], 12:[15,19], 13:[16], 14:[17], 15:[18], 16:[19]},\
                  pos = pos_dict, name="Desargues Graph")
        return G

    def FlowerSnark(self):
        """
        Returns a Flower Snark.

        A flower snark has 20 vertices. It is part of the class of
        biconnected cubic graphs with edge chromatic number = 4, known as
        snarks. (i.e.: the Petersen graph). All snarks are not Hamiltonian,
        non-planar and have Petersen graph graph minors.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algoirithm. By convention, the nodes are
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
            pos_dict[i] = [x,y]
        for i in range(20)[15:]:
            x = float(cos((pi/2) + ((2*pi)/5)*i))
            y = float(sin((pi/2) + ((2*pi)/5)*i))
            pos_dict[i] = [x,y]
        return graph.Graph({0:[1,14,15],1:[2,11],2:[3,7],3:[2,4,16],4:[5,14], \
                            5:[6,10],6:[5,7,17],8:[7,9,13],9:[10,18],11:[10,12], \
                            12:[13,19],13:[14],15:[19],16:[15,17],18:[17,19]}, \
                            pos=pos_dict, name="Flower Snark")

    def FruchtGraph(self):
        """
        Returns a Frucht Graph.

        A Frucht graph has 12 nodes and 18 edges. It is the smallest cubic
        identity graph. It is planar and it is Hamiltonian.

        This constructor is dependant on Networkx's numeric labeling.

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
            pos_dict[i] = [x,y]
        pos_dict[7] = [0,1]
        pos_dict[8] = [-1,0]
        pos_dict[9] = [0,-1]
        pos_dict[10] = [1,0]
        pos_dict[11] = [0,0]
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

        This constructor is dependant on Networkx's numeric labeling.

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
            pos_dict[i] = [x,y]
        import networkx
        G = networkx.heawood_graph()
        return graph.Graph(G, pos=pos_dict, name="Heawood graph", implementation='networkx')

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
        'p40':['p42'], 'p42':['p44'], 'p44':['p41'], 'p41':['p43'], 'p43':['p40']}, implementation='networkx' )
        for j in range(5):
            for i in range(5):
                for k in range(5):
                    con = (i+j*k)%5
                    H.add_edge(('q%d%d'%(k,con),'p%d%d'%(j,i)))
        H.name('Hoffman-Singleton graph')
        from sage.combinat.combinat import permutations
        from random import randint
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
            pos_dict[map[D[i]]] = [x,y]
        H.set_pos(pos_dict)
        return H

    def MoebiusKantorGraph(self):
        """
        Returns a Moebius-Kantor Graph.

        A Moebius-Kantor graph is a cubic symmetric graph. (See also the
        Heawood graph). It has 16 nodes and 24 edges. It is nonplanar and
        Hamiltonian. It has diameter = 4, girth = 6, and chromatic number =
        2. It is identical to the Generalized Petersen graph, P[8,3].

        PLOTTING: Upon construction, the position dictionary is filled to
        overwrite the spring-layout algorithm. By convention, the first 8
        nodes are drawn counter-clockwise in an outer circle, with the
        remaining eight drawn likewise nested in a smaller circular
        pattern. The Moebius-Kantor graph is constructed directly below
        from a dictionary with nodes as keys and entries represented the
        nodes they are connected to. Please browse this dictionary or
        display an example to further understand the plotting convention.

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
        pos_dict = {}
        for i in range(8):
            x = float(2*(cos((pi/2) + ((pi)/4)*i)))
            y = float(2*(sin((pi/2) + ((pi)/4)*i)))
            pos_dict[i] = [x,y]
        for i in range(16)[8:]:
            x = float(cos((pi/2) + ((pi)/4)*(i)))
            y = float(sin((pi/2) + ((pi)/4)*(i)))
            pos_dict[i] = [x,y]
        return graph.Graph({0:[1,7,8],1:[2,9],2:[3,10],3:[4,11],4:[5,12], \
                            5:[6,13],6:[7,14],9:[12,14],11:[8,14],13:[8,10], \
                            15:[7,10,12]}, pos=pos_dict, name="Moebius-Kantor Graph")

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

        PLOTTING: When plotting the Petersen graph with the spring-layout
        algorithm, we see that this graph is not very symmetric and thus
        the display may not be very meaningful. Efficiency of construction
        and plotting is not an issue, as the Petersen graph only has 10
        vertices.

        Our labeling convention here is to start on the outer pentagon from
        the top, moving counterclockwise. Then the nodes on the inner star,
        starting at the top and moving counterclockwise.

        EXAMPLES: We compare below the Petersen graph with the default
        spring-layout versus a planned position dictionary of [x,y]
        tuples::

            sage: petersen_spring = Graph({0:[1,4,5], 1:[0,2,6], 2:[1,3,7], 3:[2,4,8], 4:[0,3,9], 5:[0,7,8], 6:[1,8,9], 7:[2,5,9], 8:[3,5,6], 9:[4,6,7]})
            sage: petersen_spring.show() # long time
            sage: petersen_database = graphs.PetersenGraph()
            sage: petersen_database.show() # long time
        """
        pos_dict = {}
        for i in range(5):
            x = float(cos(pi/2 + ((2*pi)/5)*i))
            y = float(sin(pi/2 + ((2*pi)/5)*i))
            pos_dict[i] = [x,y]
        for i in range(10)[5:]:
            x = float(0.5*cos(pi/2 + ((2*pi)/5)*i))
            y = float(0.5*sin(pi/2 + ((2*pi)/5)*i))
            pos_dict[i] = [x,y]
        P = graph.Graph({0:[1,4,5], 1:[0,2,6], 2:[1,3,7], 3:[2,4,8], 4:[0,3,9],\
            5:[0,7,8], 6:[1,8,9], 7:[2,5,9], 8:[3,5,6], 9:[4,6,7]},\
            pos=pos_dict, name="Petersen graph")
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
        pos_dict = {0:[-1,1],1:[0,1],2:[1,1],3:[-1,0],4:[0,0],5:[1,0]}
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
            pos_dict[i] = [x,y]
        G=graph.Graph(n, name="Circulant graph ("+str(adjacency)+")")
        G._pos=pos_dict
        for v in G:
            G.add_edges([[v,(v+j)%n] for j in adjacency])
            G.add_edges([[v,(v-j)%n] for j in adjacency])
        return G



    def CompleteGraph(self, n):
        """
        Returns a complete graph on n nodes.

        A Complete Graph is a graph in which all nodes are connected to all
        other nodes.

        This constructor is dependant on vertices numbered 0 through n-1 in
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
            pos_dict[i] = [x,y]
        import networkx
        G = networkx.complete_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Complete graph", implementation='networkx')

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
            pos_dict[i] = [x,y]
        for i in range(n1+n2)[n1:]:
            x = c2*(i-n1) + c4
            y = 0
            pos_dict[i] = [x,y]
        import networkx
        import sage.graphs.bipartite_graph as bipartite_graph
        G = networkx.complete_bipartite_graph(n1,n2)
        return bipartite_graph.BipartiteGraph(G, pos=pos_dict, name="Complete bipartite graph")

    def CubeGraph(self, n):
        """
        AUTHORS:

        - Robert Miller

        PLOTTING: See commented source code.

        EXAMPLES: Plot several n-cubes in a Sage Graphics Array

        ::

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

        Use the plot options to display larger n-cubes

        ::

            sage: g = graphs.CubeGraph(9)
            sage: g.show(figsize=[12,12],vertex_labels=False, vertex_size=20) # long time
        """
        from sage.rings.integer import Integer
        # generate vertex labels:
        # n positions, 0 or 1 for each
        l = []
        for i in range(2**n):
            l.append(Integer(i).binary())
        for i in range(len(l)):
            l[i] = '0'*(n - len(l[i])) + l[i]

        # determine adjacencies:
        # adjacent vertices differ in
        # exactly one position
        d = {}
        for i in range(len(l)):
            a = []
            for j in range(n):
                if l[i][j] == '0':
                    k = '1'
                else: k = '0'
                a.append(l[i][0:j] + k + l[i][j+1:n])
            d[l[i]] = a

        # get basis vectors for projection RR^n -> RR^2
        ll = {}
        theta = float(pi/n)
        for i in range(n):
            ll[i] = (float(cos(i*theta)),float(sin(i*theta)))

        # calculate positions
        pos = {}
        for vertex in d.iterkeys():
            x = 0
            y = 0
            for i in range(n):
                x += int(vertex[i])*ll[i][0]
                y += int(vertex[i])*ll[i][1]
            pos[vertex] = [x,y]

        return graph.Graph(data=d, pos=pos, name="%d-Cube"%n, implementation='networkx')

    def BalancedTree(self, r, h):
        r"""
        Returns the perfectly balanced tree of height `h \geq 1`,
        whose root has degree `r \geq 2`.

        The number of vertices of this graph is
        `1 + r + r^2 + \cdots + r^h`, that is,
        `\frac{r^{h+1} - 1}{r - 1}`. The number of edges is one
        less than the number of vertices.

        EXAMPLE: Plot a balanced tree of height 4 with r = 3

        ::

            sage: G = graphs.BalancedTree(3, 5)
            sage: G.show()   # long time
        """
        import networkx
        return graph.Graph(networkx.balanced_tree(r, h), name="Balanced Tree")

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

        - [2] Grunbaum, B.  Convex Polytop es. New York: Wiley,
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
            G = networkx.fast_gnp_random_graph(n, p, seed)
        else:
            G = networkx.gnp_random_graph(n, p, seed)
        return graph.Graph(G)

    def RandomBarabasiAlbert(self, n, m, seed=None):
        u"""
        Return a random graph created using the Barabasi-Albert preferential
        attachment model.

        A graph with m vertices and no edges is initialized, and a graph of n
        vertices is grown by attaching new veritces each with m edges that are
        attached to existing vertices, preferentially with high degree.

        INPUT:

        - ``n`` - number of vertices in the graph

        - ``m`` - number of edges to attach from each new node

        - ``seed`` - for random number generator

        EXAMPLES:

        We show the edge list of a random graph on 6 nodes with m = 2.

        ::

            sage: graphs.RandomBarabasiAlbert(6,2).edges(labels=False)
            [(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (2, 4), (2, 5), (3, 5)]

        We plot a random graph on 12 nodes with m = 3.

        ::

            sage: ba = graphs.RandomBarabasiAlbert(12,3)
            sage: ba.show()  # long time

        We view many random graphs using a graphics array::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
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
        return graph.Graph(networkx.barabasi_albert_graph(n,m,seed))

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
            return graph.Graph(networkx.dense_gnm_random_graph(n, m, seed))
        else:
            return graph.Graph(networkx.gnm_random_graph(n, m, seed))

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
        return graph.Graph(networkx.newman_watts_strogatz_graph(n, k, p, seed))

    def RandomHolmeKim(self, n, m, p, seed=None):
        """
        Returns a random graph generated by the Holme and Kim algorithm for
        graphs with powerlaw degree distribution and approximate average
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
        return graph.Graph(networkx.powerlaw_cluster_graph(n, m, p, seed))

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
        return graph.Graph(networkx.random_lobster(n, p, q, seed))

    def RandomTreePowerlaw(self, n, gamma=3, tries=100, seed=None):
        """
        Returns a tree with a powerlaw degree distribution. Returns False
        on failure.

        From the NetworkX documentation: A trial powerlaw degree sequence
        is chosen and then elements are swapped with new elements from a
        powerlaw distribution until the sequence makes a tree (size = order
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
            return graph.Graph(networkx.random_powerlaw_tree(n, gamma, seed, tries))
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

            sage: graphs.RandomRegular(3, 8)
            Graph on 0 vertices
            sage: graphs.RandomRegular(3, 8)
            Graph on 0 vertices
            sage: graphs.RandomRegular(3, 8).edges(labels=False)
            [(0, 1), (0, 4), (0, 5), (1, 6), (1, 7), (2, 3), (2, 4), (2, 7), (3, 4), (3, 5), (5, 6), (6, 7)]

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
            return graph.Graph(networkx.random_regular_graph(d, n, seed), sparse=True)
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
        return graph.Graph(networkx.random_shell_graph(constructor, seed))

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
        return graph.Graph(networkx.configuration_model([int(i) for i in deg_sequence], seed), loops=True, multiedges=True, implementation='networkx')

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
        return graph.Graph(networkx.expected_degree_graph([int(i) for i in deg_sequence], seed), loops=True)

################################################################################
#   Graph Iterators
################################################################################

    def __call__(self, vertices, property=lambda x: True, augment='edges',
        size=None, deg_seq=None, loops=False, implementation='networkx',
        sparse=True):
        """
        Accesses the generator of isomorphism class representatives.
        Iterates over distinct, exhaustive representatives.

        INPUT:


        -  ``vertices`` - natural number

        -  ``property`` - any property to be tested on graphs
           before generation. (Ignored if deg_seq is specified.)

        -  ``augment`` - choices:

        -  ``'vertices'`` - augments by adding a vertex, and
           edges incident to that vertex. In this case, all graphs on up to
           n=vertices are generated. If for any graph G satisfying the
           property, every subgraph, obtained from G by deleting one vertex
           and only edges incident to that vertex, satisfies the property,
           then this will generate all graphs with that property. If this does
           not hold, then all the graphs generated will satisfy the property,
           but there will be some missing.

        -  ``'edges'`` - augments a fixed number of vertices by
           adding one edge In this case, all graphs on exactly n=vertices are
           generated. If for any graph G satisfying the property, every
           subgraph, obtained from G by deleting one edge but not the vertices
           incident to that edge, satisfies the property, then this will
           generate all graphs with that property. If this does not hold, then
           all the graphs generated will satisfy the property, but there will
           be some missing.

        -  ``deg_seq`` - a sequence of non-negative integers,
           or None. If specified, the generated graphs will have these
           integers for degrees. In this case property and size are both
           ignored.

        -  ``loops`` - whether to allow loops in the graph or
           not.

        -  ``implementation`` - which underlying implementation to use (see Graph?)

        -  ``sparse`` - ignored if implementation is not ``c_graph``

        EXAMPLES: Print graphs on 3 or less vertices.

        ::

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

        For more examples, see the class level documentation, or type

        ::

            sage: graphs? # not tested

        REFERENCE:

        - Brendan D. McKay, Isomorph-Free Exhaustive generation.
          Journal of Algorithms Volume 26, Issue 2, February 1998,
          pages 306-324.
        """
        from sage.graphs.graph import Graph
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

    def trees(self, vertices, augment='edges'):
        """
        Accesses the generator of trees (graphs without cycles). Iterates
        over distinct, exhaustive representatives.

        INPUT:


        -  ``vertices`` - natural number

        -  ``augment`` - choices:

        -  ``'vertices'`` - augments by adding a vertex, and
           edges incident to that vertex. In this case, all trees on up to
           n=vertices are generated.

        -  ``'edges'`` - augments a fixed number of vertices by
           adding one edge In this case, all trees on exactly n=vertices are
           generated.


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
        is_forest = lambda g: g.is_forest()
        for g in self(vertices=vertices, property=is_forest, augment=augment):
            if g.is_connected():
                yield g

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



class DiGraphGenerators():
    r"""
    A class consisting of constructors for several common digraphs,
    including orderly generation of isomorphism class representatives.

    A list of all graphs and graph structures in this database is
    available via tab completion. Type "digraphs." and then hit tab to
    see which graphs are available.

    The docstrings include educational information about each named
    digraph with the hopes that this class can be used as a reference.

    The constructors currently in this class include::

                Random Directed Graphs:
                    - RandomDirectedGN
                    - RandomDirectedGNC
                    - RandomDirectedGNR


    ORDERLY GENERATION: digraphs(vertices, property=lambda x: True,
    augment='edges', size=None)

    Accesses the generator of isomorphism class representatives.
    Iterates over distinct, exhaustive representatives.

    INPUT:


    -  ``vertices`` - natural number

    -  ``property`` - any property to be tested on digraphs
       before generation.

    -  ``augment`` - choices:

    -  ``'vertices'`` - augments by adding a vertex, and
       edges incident to that vertex. In this case, all digraphs on up to
       n=vertices are generated. If for any digraph G satisfying the
       property, every subgraph, obtained from G by deleting one vertex
       and only edges incident to that vertex, satisfies the property,
       then this will generate all digraphs with that property. If this
       does not hold, then all the digraphs generated will satisfy the
       property, but there will be some missing.

    -  ``'edges'`` - augments a fixed number of vertices by
       adding one edge In this case, all digraphs on exactly n=vertices
       are generated. If for any graph G satisfying the property, every
       subgraph, obtained from G by deleting one edge but not the vertices
       incident to that edge, satisfies the property, then this will
       generate all digraphs with that property. If this does not hold,
       then all the digraphs generated will satisfy the property, but
       there will be some missing.

    -  ``implementation`` - which underlying implementation to use (see DiGraph?)

    -  ``sparse`` - ignored if implementation is not ``c_graph``

    EXAMPLES: Print digraphs on 2 or less vertices.

    ::

        sage: for D in digraphs(2, augment='vertices'):
        ...    print D
        ...
        Digraph on 0 vertices
        Digraph on 1 vertex
        Digraph on 2 vertices
        Digraph on 2 vertices
        Digraph on 2 vertices

    Note that we can also get digraphs with underlying Cython implementation::

        sage: for D in digraphs(2, augment='vertices', implementation='c_graph'):
        ...    print D
        ...
        Digraph on 0 vertices
        Digraph on 1 vertex
        Digraph on 2 vertices
        Digraph on 2 vertices
        Digraph on 2 vertices

    Print digraphs on 3 vertices.

    ::

        sage: for D in digraphs(3):
        ...    print D
        Digraph on 3 vertices
        Digraph on 3 vertices
        ...
        Digraph on 3 vertices
        Digraph on 3 vertices

    Generate all digraphs with 4 vertices and 3 edges.

    ::

        sage: L = digraphs(4, size=3)
        sage: len(list(L))
        13

    Generate all digraphs with 4 vertices and up to 3 edges.

    ::

        sage: L = list(digraphs(4, lambda G: G.size() <= 3))
        sage: len(L)
        20
        sage: graphs_list.show_graphs(L)  # long time

    Generate all digraphs with degree at most 2, up to 5 vertices.

    ::

        sage: property = lambda G: ( max([G.degree(v) for v in G] + [0]) <= 2 )
        sage: L = list(digraphs(5, property, augment='vertices'))
        sage: len(L)
        75

    Generate digraphs on the fly: (see
    http://www.research.att.com/~njas/sequences/A000273)

    ::

        sage: for i in range(0, 5):
        ...    print len(list(digraphs(i)))
        1
        1
        3
        16
        218

    REFERENCE:

    - Brendan D. McKay, Isomorph-Free Exhaustive generation.  Journal
      of Algorithms Volume 26, Issue 2, February 1998, pages 306-324.
    """

    def ButterflyGraph(self, n, vertices='strings'):
        """
        Returns a n-dimensional butterfly graph. The vertices consist of
        pairs (v,i), where v is an n-dimensional tuple (vector) with binary
        entries (or a string representation of such) and i is an integer in
        [0..n]. A directed edge goes from (v,i) to (w,i+1) if v and w are
        identical except for possibly v[i] != w[i].

        A butterfly graph has `(2^n)(n+1)` vertices and
        `n2^{n+1}` edges.

        INPUT:


        -  ``vertices`` - 'strings' (default) or 'vectors',
           specifying whether the vertices are zero-one strings or actually
           tuples over GF(2).


        EXAMPLES::

            sage: digraphs.ButterflyGraph(2).edges(labels=False)
            [(('00', 0), ('00', 1)),
            (('00', 0), ('10', 1)),
            (('00', 1), ('00', 2)),
            (('00', 1), ('01', 2)),
            (('01', 0), ('01', 1)),
            (('01', 0), ('11', 1)),
            (('01', 1), ('00', 2)),
            (('01', 1), ('01', 2)),
            (('10', 0), ('00', 1)),
            (('10', 0), ('10', 1)),
            (('10', 1), ('10', 2)),
            (('10', 1), ('11', 2)),
            (('11', 0), ('01', 1)),
            (('11', 0), ('11', 1)),
            (('11', 1), ('10', 2)),
            (('11', 1), ('11', 2))]
            sage: digraphs.ButterflyGraph(2,vertices='vectors').edges(labels=False)
            [(((0, 0), 0), ((0, 0), 1)),
            (((0, 0), 0), ((1, 0), 1)),
            (((0, 0), 1), ((0, 0), 2)),
            (((0, 0), 1), ((0, 1), 2)),
            (((0, 1), 0), ((0, 1), 1)),
            (((0, 1), 0), ((1, 1), 1)),
            (((0, 1), 1), ((0, 0), 2)),
            (((0, 1), 1), ((0, 1), 2)),
            (((1, 0), 0), ((0, 0), 1)),
            (((1, 0), 0), ((1, 0), 1)),
            (((1, 0), 1), ((1, 0), 2)),
            (((1, 0), 1), ((1, 1), 2)),
            (((1, 1), 0), ((0, 1), 1)),
            (((1, 1), 0), ((1, 1), 1)),
            (((1, 1), 1), ((1, 0), 2)),
            (((1, 1), 1), ((1, 1), 2))]
        """
        # We could switch to Sage integers to handle arbitrary n.
        if vertices=='strings':
            if n>=31:
                raise NotImplementedError, "vertices='strings' is only valid for n<=30."
            from sage.graphs.graph_fast import binary
            butterfly = {}
            for v in xrange(2**n):
                for i in range(n):
                    w = v
                    w ^= (1 << i)   # push 1 to the left by i and xor with w
                    bv = binary(v)
                    bw = binary(w)
                    # pad and reverse the strings
                    padded_bv = ('0'*(n-len(bv))+bv)[::-1]
                    padded_bw = ('0'*(n-len(bw))+bw)[::-1]
                    butterfly[(padded_bv,i)]=[(padded_bv,i+1), (padded_bw,i+1)]
        elif vertices=='vectors':
            from sage.modules.free_module import VectorSpace
            from sage.rings.finite_field import FiniteField
            from copy import copy
            butterfly = {}
            for v in VectorSpace(FiniteField(2),n):
                for i in xrange(n):
                    w=copy(v)
                    w[i] += 1 # Flip the ith bit
                    # We must call tuple since vectors are mutable.  To obtain
                    # a vector from the tuple t, just call vector(t).
                    butterfly[(tuple(v),i)]=[(tuple(v),i+1), (tuple(w),i+1)]
        else:
            raise NotImplementedError, "vertices must be 'strings' or 'vectors'."
        return graph.DiGraph(butterfly, implementation='networkx')

    def RandomDirectedGN(self, n, kernel=lambda x:x, seed=None):
        """
        Returns a random GN (growing network) digraph with n vertices.

        The digraph is constructed by adding vertices with a link to one
        previously added vertex. The vertex to link to is chosen with a
        preferential attachment model, i.e. probability is proportional to
        degree. The default attachment kernel is a linear function of
        degree. The digraph is always a tree, so in particular it is a
        directed acyclic graph.

        INPUT:


        -  ``n`` - number of vertices.

        -  ``kernel`` - the attachment kernel

        -  ``seed`` - for the random number generator


        EXAMPLE::

            sage: D = digraphs.RandomDirectedGN(25)
            sage: D.edges(labels=False)
            [(1, 0), (2, 0), (3, 1), (4, 0), (5, 0), (6, 1), (7, 0), (8, 3), (9, 0), (10, 8), (11, 3), (12, 9), (13, 8), (14, 0), (15, 11), (16, 11), (17, 5), (18, 11), (19, 6), (20, 5), (21, 14), (22, 5), (23, 18), (24, 11)]
            sage: D.show()  # long time

        REFERENCE:

        - [1] Krapivsky, P.L. and Redner, S. Organization of Growing
          Random Networks, Phys. Rev. E vol. 63 (2001), p. 066123.
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return graph.DiGraph(networkx.gn_graph(n, kernel, seed))

    def RandomDirectedGNC(self, n, seed=None):
        """
        Returns a random GNC (growing network with copying) digraph with n
        vertices.

        The digraph is constructed by adding vertices with a link to one
        previously added vertex. The vertex to link to is chosen with a
        preferential attachment model, i.e. probability is proportional to
        degree. The new vertex is also linked to all of the previously
        added vertex's successors.

        INPUT:


        -  ``n`` - number of vertices.

        -  ``seed`` - for the random number generator


        EXAMPLE::

            sage: D = digraphs.RandomDirectedGNC(25)
            sage: D.edges(labels=False)
            [(1, 0), (2, 0), (2, 1), (3, 0), (4, 0), (4, 1), (5, 0), (5, 1), (5, 2), (6, 0), (6, 1), (7, 0), (7, 1), (7, 4), (8, 0), (9, 0), (9, 8), (10, 0), (10, 1), (10, 2), (10, 5), (11, 0), (11, 8), (11, 9), (12, 0), (12, 8), (12, 9), (13, 0), (13, 1), (14, 0), (14, 8), (14, 9), (14, 12), (15, 0), (15, 8), (15, 9), (15, 12), (16, 0), (16, 1), (16, 4), (16, 7), (17, 0), (17, 8), (17, 9), (17, 12), (18, 0), (18, 8), (19, 0), (19, 1), (19, 4), (19, 7), (20, 0), (20, 1), (20, 4), (20, 7), (20, 16), (21, 0), (21, 8), (22, 0), (22, 1), (22, 4), (22, 7), (22, 19), (23, 0), (23, 8), (23, 9), (23, 12), (23, 14), (24, 0), (24, 8), (24, 9), (24, 12), (24, 15)]
            sage: D.show()  # long time

        REFERENCE:

        - [1] Krapivsky, P.L. and Redner, S. Network Growth by
          Copying, Phys. Rev. E vol. 71 (2005), p. 036118.
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return graph.DiGraph(networkx.gnc_graph(n, seed))

    def RandomDirectedGNP(self, n, p):
        r"""
        Returns a random digraph on `n` nodes. Each edge is
        inserted independently with probability `p`.

        REFERENCES:

        - [1] P. Erdos and A. Renyi, On Random Graphs, Publ.  Math. 6,
          290 (1959).

        - [2] E. N. Gilbert, Random Graphs, Ann. Math.  Stat., 30,
          1141 (1959).

        PLOTTING: When plotting, this graph will use the default
        spring-layout algorithm, unless a position dictionary is
        specified.

        EXAMPLE::

            sage: digraphs.RandomDirectedGNP(10, .2).num_verts()
            10
        """
        from random import random
        D = graph.DiGraph(n)
        for i in xrange(n):
            for j in xrange(i):
                if random() < p:
                    D.add_edge(i,j)
            for j in xrange(i+1,n):
                if random() < p:
                    D.add_edge(i,j)
        return D

    def RandomDirectedGNR(self, n, p, seed=None):
        """
        Returns a random GNR (growing network with redirection) digraph
        with n vertices and redirection probability p.

        The digraph is constructed by adding vertices with a link to one
        previously added vertex. The vertex to link to is chosen uniformly.
        With probability p, the arc is instead redirected to the successor
        vertex. The digraph is always a tree.

        INPUT:


        -  ``n`` - number of vertices.

        -  ``p`` - redirection probability

        -  ``seed`` - for the random number generator.


        EXAMPLE::

            sage: D = digraphs.RandomDirectedGNR(25, .2)
            sage: D.edges(labels=False)
            [(1, 0), (2, 0), (2, 1), (3, 0), (4, 0), (4, 1), (5, 0), (5, 1), (5, 2), (6, 0), (6, 1), (7, 0), (7, 1), (7, 4), (8, 0), (9, 0), (9, 8), (10, 0), (10, 1), (10, 2), (10, 5), (11, 0), (11, 8), (11, 9), (12, 0), (12, 8), (12, 9), (13, 0), (13, 1), (14, 0), (14, 8), (14, 9), (14, 12), (15, 0), (15, 8), (15, 9), (15, 12), (16, 0), (16, 1), (16, 4), (16, 7), (17, 0), (17, 8), (17, 9), (17, 12), (18, 0), (18, 8), (19, 0), (19, 1), (19, 4), (19, 7), (20, 0), (20, 1), (20, 4), (20, 7), (20, 16), (21, 0), (21, 8), (22, 0), (22, 1), (22, 4), (22, 7), (22, 19), (23, 0), (23, 8), (23, 9), (23, 12), (23, 14), (24, 0), (24, 8), (24, 9), (24, 12), (24, 15)]
            sage: D.show()  # long time

        REFERENCE:

        - [1] Krapivsky, P.L. and Redner, S. Organization of Growing
          Random Networks, Phys. Rev. E vol. 63 (2001), p. 066123.
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return graph.DiGraph(networkx.gnc_graph(n, seed))

################################################################################
#   DiGraph Iterators
################################################################################

    def __call__(self, vertices, property=lambda x: True, augment='edges', size=None, implementation='networkx', sparse=True):
        """
        Accesses the generator of isomorphism class representatives.
        Iterates over distinct, exhaustive representatives.

        INPUT:


        -  ``vertices`` - natural number

        -  ``property`` - any property to be tested on digraphs
           before generation.

        -  ``augment`` - choices:

        -  ``'vertices'`` - augments by adding a vertex, and
           edges incident to that vertex. In this case, all digraphs on up to
           n=vertices are generated. If for any digraph G satisfying the
           property, every subgraph, obtained from G by deleting one vertex
           and only edges incident to that vertex, satisfies the property,
           then this will generate all digraphs with that property. If this
           does not hold, then all the digraphs generated will satisfy the
           property, but there will be some missing.

        -  ``'edges'`` - augments a fixed number of vertices by
           adding one edge In this case, all digraphs on exactly n=vertices
           are generated. If for any graph G satisfying the property, every
           subgraph, obtained from G by deleting one edge but not the vertices
           incident to that edge, satisfies the property, then this will
           generate all digraphs with that property. If this does not hold,
           then all the digraphs generated will satisfy the property, but
           there will be some missing.

        -  ``implementation`` - which underlying implementation to use (see DiGraph?)

        -  ``sparse`` - ignored if implementation is not ``c_graph``

        EXAMPLES: Print digraphs on 2 or less vertices.

        ::

            sage: for D in digraphs(2, augment='vertices'):
            ...    print D
            ...
            Digraph on 0 vertices
            Digraph on 1 vertex
            Digraph on 2 vertices
            Digraph on 2 vertices
            Digraph on 2 vertices

        Print digraphs on 3 vertices.

        ::

            sage: for D in digraphs(3):
            ...    print D
            Digraph on 3 vertices
            Digraph on 3 vertices
            ...
            Digraph on 3 vertices
            Digraph on 3 vertices

        For more examples, see the class level documentation, or type

        ::

            sage: digraphs? # not tested

        REFERENCE:

        - Brendan D. McKay, Isomorph-Free Exhaustive generation.
          Journal of Algorithms Volume 26, Issue 2, February 1998,
          pages 306-324.
        """
        from sage.graphs.graph import DiGraph
        if size is not None:
            extra_property = lambda x: x.size() == size
        else:
            extra_property = lambda x: True
        if augment == 'vertices':
            g = DiGraph(implementation=implementation, sparse=sparse)
            for gg in canaug_traverse_vert(g, [], vertices, property, dig=True, implementation=implementation, sparse=sparse):
                if extra_property(gg):
                    yield gg
        elif augment == 'edges':
            g = DiGraph(vertices, implementation=implementation, sparse=sparse)
            gens = []
            for i in range(vertices-1):
                gen = range(i)
                gen.append(i+1); gen.append(i)
                gen += range(i+2, vertices)
                gens.append(gen)
            for gg in canaug_traverse_edge(g, gens, property, dig=True, implementation=implementation, sparse=sparse):
                if extra_property(gg):
                    yield gg
        else:
            raise NotImplementedError()

def canaug_traverse_vert(g, aut_gens, max_verts, property, dig=False, loops=False, implementation='networkx', sparse=True):
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
    from sage.graphs.graph_fast import binary
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

def canaug_traverse_edge(g, aut_gens, property, dig=False, loops=False, implementation='networkx', sparse=True):
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
    from sage.graphs.graph_fast import binary
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

            m_z = z.copy()
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
digraphs = DiGraphGenerators()




