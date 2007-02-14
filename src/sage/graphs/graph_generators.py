r"""
A collection of constructors of common graphs.

USE:

    To see a list of all graph constructors, type "graphs." and then press the
    tab key.  The documentation for each constructor includes information about
    each graph, which provides a useful reference.

PLOTTING:
    All graphs (i.e., networks) have an associated SAGE graphics object,
    which you can display:

        sage: G = graphs.WheelGraph(15)
        sage: P = G.plot()
        sage.: P.show()

    If you create a graph in SAGE using the \code{Graph} command, then
    plot that graph, the positioning of nodes is determined using the
    spring-layout algorithm.  For the special graph constructors,
    which you get using \code{graphs.[tab]}, the positions are preset.
    For example, consider the Petersen graph with default node
    positioning vs. the Petersen graph constructed by this database:

        sage: petersen_spring = Graph({0:[1,4,5], 1:[0,2,6], 2:[1,3,7], 3:[2,4,8], 4:[0,3,9], 5:[0,7,8], 6:[1,8,9], 7:[2,5,9], 8:[3,5,6], 9:[4,6,7]})
        sage.: petersen_spring.show()
        sage: petersen_database = graphs.PetersenGraph()
        sage.: petersen_database.show()

    For all the constructors in this database (except the octahedral,
    dodecahedral, random and empty graphs), the position dictionary
    is filled in, instead of using the spring-layout algorithm.

    For further visual examples and explanation, see the docstrings
    below, particularly for CycleGraph, StarGraph, WheelGraph,
    CompleteGraph and CompleteBipartiteGraph.

ORGANIZATION:
    The constructors available in this database are organized as follows:
    \begin{verbatim}
        Basic Structures:
            - BarbellGraph
            - BullGraph
            - CircularLadderGraph
            - CycleGraph
            - DiamondGraph
            - DodecahedralGraph
            - EmptyGraph
            - Grid2dGraph
            - HouseGraph
            - HouseXGraph
            - KrackhardtKiteGraph
            - LadderGraph
            - LollipopGraph
            - OctahedralGraph
            - PathGraph
            - StarGraph
            - TetrahedralGraph
            - WheelGraph
        Named Graphs:
            - PetersenGraph
        Families of Graphs:
            - CompleteGraph
            - CompleteBipartiteGraph
            - CubeGraph
            - RandomGNP
            - RandomGNPFast
    \end{verbatim}

AUTHORS:
    -- Robert Miller (2006-11-05): initial version - empty, random,
       petersen
    -- Emily Kirkman (2006-11-12): basic structures, node positioning for
       all constructors
    -- Emily Kirkman (2006-11-19): docstrings, examples
    -- William Stein (2006-12-05): Editing.
    -- Robert Miller (2007-01-16): Cube generation and plotting
    -- Emily Kirkman (2007-01-16): more basic structures, docstrings
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

class GraphGenerators():
    r"""
    A class consisting of constructors for several common graphs.

    A list of all graphs and graph structures in this database is available
    via tab completion. Type "graphs." and then hit tab to see which graphs
    are available.

    The docstrings include educational information about each named graph
    with the hopes that this class can be used as a reference.

    For all the constructors in this class (except the octahedral,
    dodecahedral, random and empty graphs), the position dictionary
    is filled to override the spring-layout algorithm.

    The constructors currently in this class include:
    \begin{verbatim}
        Basic Structures:
             - BarbellGraph
             - BullGraph
             - CircularLadderGraph
             - CycleGraph
             - DiamondGraph
             - DodecahedralGraph
             - EmptyGraph
             - Grid2dGraph
             - HouseGraph
             - HouseXGraph
             - KrackhardtKiteGraph
             - LadderGraph
             - LollipopGraph
             - OctahedralGraph
             - PathGraph
             - StarGraph
             - TetrahedralGraph
             - WheelGraph
        Named Graphs:
            - PetersenGraph
        Families of Graphs:
            - CompleteGraph
            - CompleteBipartiteGraph
            - CubeGraph
            - RandomGNP
            - RandomGNPFast
    \end{verbatim}
    """

################################################################################
#   Basic Structures
################################################################################

    def BarbellGraph(self, n1, n2):
        """
        Returns a barbell graph with 2*n1 + n2 nodes.
        n1 must be greater than or equal to 2.

        A barbell graph is a basic structure that consists of a path graph of order
        n2 connecting two complete graphs of order n1 each.

        This constructor depends on NetworkX numeric labels.  In this case, the
        (n1)th node connects to the path graph from one complete graph and the
        (n1+n2+1)th node connects to the path graph from the other complete graph.

        PLOTTING:
        Upon construction, the position dictionary is filled to override
        the spring-layout algorithm. By convention, each barbell graph will
        be displayed with the two complete graphs in the lower-left and
        upper-right corners, with the path graph connecting diagonally
        between the two.  Thus the (n1)th node will be drawn at a 45 degree
        angle from the horizontal right center of the first complete graph,
        and the (n1+n2+1)th node will be drawn 45 degrees below the left
        horizontal center of the second complete graph.

        EXAMPLES:
            # Construct and show a barbell graph
            # Bar = 4, Bells = 9
            sage: g = graphs.BarbellGraph(9,4)
            sage.: g.show()

            # Create several barbell graphs in a SAGE graphics array
            sage: g = []
            sage: j = []
            sage: for i in range(6):
            ...    k = graphs.BarbellGraph(i+2,4)
            ...    g.append(k)
            ...
            sage: for i in range(2):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(node_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()
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

        A bull graph is named for its shape.  It's a triangle
        with horns.

        This constructor depends on NetworkX numeric labeling.

        PLOTTING:
        Upon construction, the position dictionary is filled to override
        the spring-layout algorithm.  By convention, the bull graph is
        drawn as a triangle with the first node (0) on the bottom.  The
        second and third nodes (1 and 2) complete the triangle.  Node 3
        is the horn connected to 1 and node 4 is the horn connected to
        node 2.

        EXAMPLES:
            # Construct and show a bull graph
            sage: g = graphs.BullGraph()
            sage.: g.show()
        """
        pos_dict = {0:[0,0],1:[-1,1],2:[1,1],3:[-2,2],4:[2,2]}
        import networkx
        G = networkx.bull_graph()
        return graph.Graph(G, pos=pos_dict, name="Bull Graph")


    def CircularLadderGraph(self, n):
        """
        Returns a circular ladder graph with 2*n nodes.

        A Circular ladder graph is a ladder graph that is connected at the
        ends, i.e.: a ladder bent around so that top meets bottom.  Thus it
        can be described as two parrallel cycle graphs connected at each
        corresponding node pair.

        This constructor depends on NetworkX numeric labels.

        PLOTTING:
        Upon construction, the position dictionary is filled to override the
        spring-layout algorithm.  By convention, the circular ladder graph is
        displayed as an inner and outer cycle pair, with the first n nodes
        drawn on the inner circle.  The first (0) node is drawn at the top
        of the inner-circle, moving clockwise after that.  The outer circle
        is drawn with the (n+1)th node at the top, then counterclockwise as
        well.

        EXAMPLES:
            # Construct and show a circular ladder graph with 26 nodes
            sage: g = graphs.CircularLadderGraph(13)
            sage.: g.show()

            # Create several circular ladder graphs in a SAGE graphics array
            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.CircularLadderGraph(i+3)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(node_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()
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

    def CycleGraph(self, n):
        r"""
        Returns a cycle graph with n nodes.

        A cycle graph is a basic structure which is also typically called
        an n-gon.

        This constructor is dependant on vertices numbered 0 through n-1
        in NetworkX \code{cycle_graph()}

        PLOTTING:
        Upon construction, the position dictionary is filled to override
        the spring-layout algorithm. By convention, each cycle graph will
        be displayed with the first (0) node at the top, with the rest
        following in a counterclockwise manner.

        The cycle graph is a good opportunity to compare efficiency of
        filling a position dictionary vs. using the spring-layout algorithm
        for plotting.  Because the cycle graph is very symmetric, the
        resulting plots should be similar (in cases of small n).

        Filling the position dictionary in advance adds O(n) to the
        constructor.  Feel free to race the constructors below in the
        examples section.  The much larger difference is the time added
        by the spring-layout algorithm when plotting.  (Also shown in the
        example below).  The spring  model is typically described as O(n^3),
        as appears to be the case in the NetworkX source code.

        EXAMPLES:

            # Compare the constructors (results will vary)
            sage: import networkx
            sage.: time n = networkx.cycle_graph(3989); spring3989 = Graph(n)
            # CPU time: 0.05 s,  Wall time: 0.07 s
            sage.: time posdict3989 = graphs.CycleGraph(3989)
            # CPU time: 5.18 s,  Wall time: 6.17 s

            # Compare the plotting speeds (results will vary)
            sage: n = networkx.cycle_graph(23)
            sage: spring23 = Graph(n)
            sage: posdict23 = graphs.CycleGraph(23)
            sage.: time spring23.show()
            # CPU time: 2.04 s,  Wall time: 2.72 s
            sage.: time posdict23.show()
            # CPU time: 0.57 s,  Wall time: 0.71 s

        We next view many cycle graphs as a SAGE graphics array.
        First we use the \code{CycleGraph} constructor, which fills in
        the position dictionary:

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.CycleGraph(i+3)
            ...    g.append(k)
            ...
            sage.: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(node_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage.: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()

        Compare to plotting with the spring-layout algorithm:
            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    spr = networkx.cycle_graph(i+3)
            ...    k = Graph(spr)
            ...    g.append(k)
            ...
            sage.: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(node_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage.: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()
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

        PLOTTING:
        Upon construction, the position dictionary is filled to override
        the spring-layout algorithm.  By convention, the diamond graph
        is drawn as a diamond, with the first node on top, second on the
        left, third on the right, and fourth on the bottom; with the
        second and third node connected.

        EXAMPLES:
            # Construct and show a diamond graph
            sage: g = graphs.DiamondGraph()
            sage.: g.show()
        """
        pos_dict = {0:[0,1],1:[-1,0],2:[1,0],3:[0,-1]}
        import networkx
        G = networkx.diamond_graph()
        return graph.Graph(G, pos=pos_dict, name="Diamond Graph")

    def DodecahedralGraph(self):
        """
        Returns a Dodecahedral graph (with 20 nodes)

        The dodecahedral graph is cubic symmetric, so the spring-layout
        algorithm will be very effective for display.

        PLOTTING:
        The Dodecahedral graph should be viewed in 3 dimensions.  We
        chose to use the default spring-layout algorithm here, so that
        multiple iterations might yield a different point of reference for
        the user.  We hope to add rotatable, 3-dimensional viewing in
        the future.  In such a case, a string argument will be added to select
        the flat spring-layout over a future implementation.

        EXAMPLES:
            # Construct and show a Dodecahdedral graph
            sage: g = graphs.DodecahedralGraph()
            sage.: g.show()

            # Create several dodecahedral graphs in a SAGE graphics array
            # They will be drawn differently due to the use of the spring-layout algorithm
            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.DodecahedralGraph()
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(node_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()
        """
        import networkx
        G = networkx.dodecahedral_graph()
        return graph.Graph(G, name="Dodecahedral")

    def EmptyGraph(self):
        """
        Returns an empty graph (0 nodes and 0 edges).

        This is useful for constructing graphs by adding edges and
        vertices individually or in a loop.

        PLOTTING:
        When plotting, this graph will use the default spring-layout
        algorithm, unless a position dictionary is specified.

        EXAMPLES:
            # Add one vertex to an empty graph and then show:
            sage: empty1 = graphs.EmptyGraph()
            sage: empty1.add_vertex()
            sage.: empty1.show()

            # Use for loops to build a graph from an empty graph:
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
            sage.: empty2.show()
        """
        return graph.Graph()

    def Grid2dGraph(self, n1, n2):
        """
        Returns a 2-dimensional grid graph with n1*n2 nodes (n1 rows and n2 columns).

        A 2d grid graph resembles a 2 dimensional grid.  All inner nodes are
        connected to their 4 neighbors.  Outer (non-corner) nodes are connected
        to their 3 neighbors.  Corner nodes are connected to their 2 neighbors.

        This constructor depends on NetworkX numeric labels.

        PLOTTING:
        Upon construction, the position dictionary is filled to override the
        spring-layout algorithm.  By convention, nodes are labelled in
        (row, column) pairs with (0, 0) in the top left corner.  Edges will
        always be horizontal and vertical - another advantage of filling the
        position dictionary.

        EXAMPLES:
            # Construct and show a grid 2d graph
            # Rows = 5, Columns = 7
            sage: g = graphs.Grid2dGraph(5,7)
            sage.: g.show()
        """
        pos_dict = {}
        for i in range(n1):
            y = -i
            for j in range(n2):
                x = j
                pos_dict[i,j] = [x,y]
        import networkx
        G = networkx.grid_2d_graph(n1,n2)
        return graph.Graph(G, pos=pos_dict, name="2D Grid Graph")

    def HouseGraph(self):
        """
        Returns a house graph with 5 nodes.

        A house graph is named for its shape.  It is a triange (roof)
        over a square (walls).

        This constructor depends on NetworkX numeric labeling.

        PLOTTING:
        Upon construction, the position dictionary is filled to override
        the spring-layout algorithm.  By convention, the house graph is
        drawn with the first node in the lower-left corner of the house,
        the second in the lower-right corner of the house.  The third
        node is in the upper-left corner connecting the roof to the wall,
        and the fourth is in the upper-right corner connecting the roof
        to the walll.  The fifth node is the top of the roof, connected
        only to the third and fourth.

        EXAMPLES:
            # Construct and show a house graph
            sage: g = graphs.HouseGraph()
            sage.: g.show()
        """
        pos_dict = {0:[-1,0],1:[1,0],2:[-1,1],3:[1,1],4:[0,2]}
        import networkx
        G = networkx.house_graph()
        return graph.Graph(G, pos=pos_dict, name="House Graph")

    def HouseXGraph(self):
        """
        Returns a house X graph with 5 nodes.

        A house X graph is a house graph with two additional edges.
        The upper-right corner is connected to the lower-left.  And
        the upper-left corner is connected to the lower-right.

        This constructor depends on NetworkX numeric labeling.

        PLOTTING:
        Upon construction, the position dictionary is filled to override
        the spring-layout algorithm.  By convention, the house X graph is
        drawn with the first node in the lower-left corner of the house,
        the second in the lower-right corner of the house.  The third
        node is in the upper-left corner connecting the roof to the wall,
        and the fourth is in the upper-right corner connecting the roof
        to the walll.  The fifth node is the top of the roof, connected
        only to the third and fourth.

        EXAMPLES:
            # Construct and show a house X graph
            sage: g = graphs.HouseXGraph()
            sage.: g.show()
        """
        pos_dict = {0:[-1,0],1:[1,0],2:[-1,1],3:[1,1],4:[0,2]}
        import networkx
        G = networkx.house_x_graph()
        return graph.Graph(G, pos=pos_dict, name="House Graph")

    def KrackhardtKiteGraph(self):
        """
        Returns a Krackhardt kite graph with 10 nodes.

        The Krackhardt kite graph was originally developed by David
        Krackhardt for the purpose of studying social networks.  It
        is used to show the distinction between:  degree centrality,
        betweeness centrality, and closeness centrality.  For more
        information read the plotting section below in conjunction
        with the example.

        REFERENCES:
            Kreps, V. (2002). "Social Network Analysis". [Online] Available: http://www.fsu.edu/~spap/water/network/intro.htm [2007, January 17]

        This constructor depends on NetworkX numeric labeling.

        PLOTTING:
        Upon construction, the position dictionary is filled to override
        the spring-layout algorithm.  By convention, the graph is drawn
        left to right, in top to bottom row sequence of [2, 3, 2, 1, 1, 1]
        nodes on each row.  This places the fourth node (3) in the center of
        the kite, with the highest degree.  But the fourth node only connects
        nodes that are otherwise connected, or those in its clique (i.e.:
        Degree Centrality).  The eigth (7) node is where the kite meets the
        tail.  It has degree = 3, less than the average, but is the only
        connection between the kite and tail (i.e.: Betweenness Centrality).
        The sixth and seventh nodes (5 and 6) are drawn in the third row and
        have degree = 5.  These nodes have the shortest path to all other nodes
        in the graph (i.e.: Closeness Centrality).  Please execute the
        example for visualization.

        EXAMPLE:
            # Construct and show a Krackhardt kite graph
            sage: g = graphs.KrackhardtKiteGraph()
            sage.: g.show()
        """
        pos_dict = {0:[-1,4],1:[1,4],2:[-2,3],3:[0,3],4:[2,3],5:[-1,2],6:[1,2],7:[0,1],8:[0,0],9:[0,-1]}
        import networkx
        G = networkx.krackhardt_kite_graph()
        return graph.Graph(G, pos=pos_dict, name="Krackhardt Kite Graph")

    def LadderGraph(self, n):
        """
        Returns a ladder graph with 2*n nodes.

        A ladder graph is a basic structure that is typically displayed as a
        ladder, i.e.:  two parallel path graphs connected at each corresponding
        node pair.

        This constructor depends on NetworkX numeric labels.

        PLOTTING:
        Upon construction, the position dictionary is filled to override
        the spring-layout algorithm. By convention, each ladder graph will
        be displayed horizontally, with the first n nodes displayed left to
        right on the top horizontal line.

        EXAMPLES:
            # Construct and show a ladder graph with 14 nodes
            sage: g = graphs.LadderGraph(7)
            sage.: g.show()

            # Create several ladder graphs in a SAGE graphics array
            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.LadderGraph(i+2)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(node_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()
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
        graph (order n1).  (A barbell graph minus one of the bells).

        This constructor depends on NetworkX numeric labels.

        PLOTTING:
        Upon construction, the position dictionary is filled to override
        the spring-layout algorithm.  By convention, the complete graph
        will be drawn in the lower-left corner with the (n1)th node at a
        45 degree angle above the right horizontal center of the complete
        graph, leading directly into the path graph.

        EXAMPLES:
            # Construct and show a lollipop graph
            # Candy = 13, Stick = 4
            sage: g = graphs.LollipopGraph(13,4)
            sage.: g.show()

            # Create several lollipop graphs in a SAGE graphics array
            sage: g = []
            sage: j = []
            sage: for i in range(6):
            ...    k = graphs.LollipopGraph(i+3,4)
            ...    g.append(k)
            ...
            sage: for i in range(2):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(node_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()
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

    def OctahedralGraph(self):
        """
        Returns an Octahedral graph (with 6 nodes).

        The Octahedral graph is the line graph of the Tetrahedral.  The
        octahedral is symmetric, so the spring-layout algorithm will be
        very effective for display.

        PLOTTING:
        The Octahedral graph should be viewed in 3 dimensions.  We
        chose to use the default spring-layout algorithm here, so that
        multiple iterations might yield a different point of reference for
        the user.  We hope to add rotatable, 3-dimensional viewing in
        the future.  In such a case, a string argument will be added to select
        the flat spring-layout over a future implementation.

        EXAMPLES:
            # Construct and show an Octahedral graph
            sage: g = graphs.OctahedralGraph()
            sage.: g.show()

            # Create several octahedral graphs in a SAGE graphics array
            # They will be drawn differently due to the use of the spring-layout algorithm
            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.OctahedralGraph()
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(node_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()
        """
        import networkx
        G = networkx.octahedral_graph()
        return graph.Graph(G, name="Octahedral")

    def PathGraph(self, n, pos=None):
        """
        Returns a path graph with n nodes.
        Pos argument takes a string which is either 'circle' or 'line',
        (otherwise the default is used).  See the plotting section below
        for more detail.

        A path graph is a graph where all inner nodes are connected
        to their two neighbors and the two end-nodes are connected to
        their one inner neighbors.  (i.e.: a cycle graph without the
        first and last node connected).

        This constructor depends on NetworkX numeric labels.

        PLOTTING:
        Upon construction, the position dictionary is filled to override
        the spring-layout algorithm.  By convention, the graph may be
        drawn in one of two ways:  The 'line' argument will draw the graph
        in a horizontal line (left to right) if there are less than 11
        nodes.  Otherwise the 'line' argument will append horizontal lines
        of length 10 nodes below, alternating left to right and right to
        left.  The 'circle' argument will cause the graph to be drawn in
        a cycle-shape, with the first node at the top and then about the
        circle in a clockwise manner.  By default (without an appropriate
        string argument) the graph will be drawn as a 'circle' if
        10 < n < 41 and as a 'line' for all other n.

        EXAMPLES:
            # Show default drawing by size:
            # 'line': n < 11
            sage: p = graphs.PathGraph(10)
            sage.: p.show()

            # 'circle': 10 < n < 41
            sage: q = graphs.PathGraph(25)
            sage.: q.show()

            # 'line': n > 40
            sage: r = graphs.PathGraph(55)
            sage.: r.show()

            # Override the default drawing:
            sage: s = graphs.PathGraph(5,'circle')
            sage.: s.show()
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

        PLOTTING:
        Upon construction, the position dictionary is filled to override
        the spring-layout algorithm. By convention, each star graph will
        be displayed with the first (0) node in the center, the second
        node (1) at the top, with the rest following in a counterclockwise
        manner.  (0) is the node connected to all other nodes.

        The star graph is a good opportunity to compare efficiency of
        filling a position dictionary vs. using the spring-layout algorithm
        for plotting.  As far as display, the spring-layout should push all
        other nodes away from the (0) node, and thus look very similar to
        this constructor's positioning.

        Filling the position dictionary in advance adds O(n) to the
        constructor.  Feel free to race the constructors below in the
        examples section.  The much larger difference is the time added
        by the spring-layout algorithm when plotting.  (Also shown in the
        example below).  The spring model is typically described as O(n^3),
        as appears to be the case in the NetworkX source code.

        EXAMPLES:
            sage: import networkx

        Compare the constructors (results will vary)
            sage.: time n = networkx.star_graph(3989); spring3989 = Graph(n)
            # CPU time: 0.08 s,  Wall time: 0.10 s
            sage.: time posdict3989 = graphs.StarGraph(3989)
            # CPU time: 5.43 s,  Wall time: 7.41 s

        Compare the plotting speeds (results will vary)
            sage: n = networkx.star_graph(23)
            sage: spring23 = Graph(n)
            sage: posdict23 = graphs.StarGraph(23)
            sage.: time spring23.show()
            # CPU time: 2.31 s,  Wall time: 3.14 s
            sage.: time posdict23.show()
            # CPU time: 0.68 s,  Wall time: 0.80 s

        View many star graphs as a SAGE Graphics Array

        With this constructor (i.e., the position dictionary filled)
            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.StarGraph(i+3)
            ...    g.append(k)
            ...
            sage.: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(node_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage.: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()

        Compared to plotting with the spring-layout algorithm
            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    spr = networkx.star_graph(i+3)
            ...    k = Graph(spr)
            ...    g.append(k)
            ...
            sage.: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(node_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage.: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()
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

    def TetrahedralGraph(self):
        """
        Returns a Tetrahedral graph (with 4 nodes).

        A tetrahedron is a 4-sided triangular pyramid.  This graph is equivalent
        to a Wheel graph with 4 nodes and also a Complete graph on four nodes.
        (See examples below).

        This constructor depends on NetworkX numeric labeling

        PLOTTING:
        Upon construction, the position dictionary is filled to override
        the spring-layout algorithm.  By convention, the graph is drawn
        with the first node at the vertical top (i.e.: positive z-axis), the
        second node in the bottom left (i.e.: positive x-axis), the third node
        in the center (i.e.: origin), and the fourth node on the right (i.e.:
        positive y-axis).  The references to axes are by the right-handed
        Cartesian coordinate system and are meant only to describe the
        appearance.  The position dictionary is actually filled with only
        (x,y) pairs.

        EXAMPLES:
            # Construct and show a Tetrahedral graph
            sage: g = graphs.TetrahedralGraph()
            sage: g.save('sage.png')

            # The following example requires networkx:
            sage: import networkx as NX

            # Compare this Tetrahedral, Wheel(4), Complete(4), and the
            # Tetrahedral plotted with the spring-layout algorithm below
            # in a SAGE graphics array:
            sage: tetra_pos = graphs.TetrahedralGraph()
            sage: tetra_spring = Graph(NX.tetrahedral_graph())
            sage: wheel = graphs.WheelGraph(4)
            sage: complete = graphs.CompleteGraph(4)
            sage: g = [tetra_pos, tetra_spring, wheel, complete]
            sage: j = []
            sage: for i in range(2):
            ...    n = []
            ...    for m in range(2):
            ...        n.append(g[i + m].plot(node_size=50, vertex_labels=False))
            ...    j.append(n)
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.save('sage.png')
        """
        pos_dict = {0:[0,1],1:[-.71,-.71],2:[0,0],3:[1.3,0]}
        import networkx
        G = networkx.tetrahedral_graph()
        return graph.Graph(G, pos=pos_dict, name="Tetrahedral")

    def WheelGraph(self, n):
        """
        Returns a Wheel graph with n nodes.

        A Wheel graph is a basic structure where one node is connected to
        all other nodes and those (outer) nodes are connected cyclically.

        This constructor depends on NetworkX numeric labels.

        PLOTTING:
        Upon construction, the position dictionary is filled to override
        the spring-layout algorithm. By convention, each wheel graph will
        be displayed with the first (0) node in the center, the second node
        at the top, and the rest following in a counterclockwise manner.

        With the wheel graph, we see that it doesn't take a very large n
        at all for the spring-layout to give a counter-intuitive display.
        (See Graphics Array examples below).

        Filling the position dictionary in advance adds O(n) to the
        constructor.  Feel free to race the constructors below in the
        examples section.  The much larger difference is the time added
        by the spring-layout algorithm when plotting.  (Also shown in the
        example below).  The spring model is typically described as O(n^3),
        as appears to be the case in the NetworkX source code.

        EXAMPLES:
        We view many wheel graphs with a SAGE Graphics Array, first
        with this constructor (i.e., the position dictionary filled):
            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.WheelGraph(i+3)
            ...    g.append(k)
            ...
            sage.: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(node_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage.: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()

        Next, using the spring-layout algorithm:
            sage: import networkx
            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    spr = networkx.wheel_graph(i+3)
            ...    k = Graph(spr)
            ...    g.append(k)
            ...
            sage.: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(node_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage.: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()

        Compare the constructors (results will vary):
            sage.: time n = networkx.wheel_graph(3989); spring3989 = Graph(n)
            # CPU time: 0.07 s,  Wall time: 0.09 s
            sage.: time posdict3989 = graphs.WheelGraph(3989)
            # CPU time: 5.99 s,  Wall time: 8.74 s

            # Compare the plotting speeds (results will vary)
            sage: n = networkx.wheel_graph(23)
            sage: spring23 = Graph(n)
            sage: posdict23 = graphs.WheelGraph(23)
            sage.: time spring23.show()
            # CPU time: 2.24 s,  Wall time: 3.00 s
            sage.: time posdict23.show()
            # CPU time: 0.68 s,  Wall time: 1.14 s
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
#   Named Graphs
################################################################################

    def PetersenGraph(self):
        """
        The Petersen Graph is a named graph that consists of 10 vertices
        and 14 edges, usually drawn as a five-point star embedded in a
        pentagon.

        The Petersen Graph is a common counterexample.  For example, it is
        not Hamiltonian.

        PLOTTING:
        When plotting the Petersen graph with the spring-layout algorithm,
        we see that this graph is not very symmetric and thus the display
        may not be very meaningful. Efficiency of construction and plotting
        is not an issue, as the Petersen graph only has 10 vertices and 14
        edges.

        Our labeling convention here is to start on the outer pentagon from
        the top, moving counterclockwise. Then the nodes on the inner star,
        starting at the top and moving counterclockwise.

        EXAMPLES:
        We compare below the Petersen graph with the default spring-layout
        versus a planned position dictionary of [x,y] tuples:
            sage: petersen_spring = Graph({0:[1,4,5], 1:[0,2,6], 2:[1,3,7], 3:[2,4,8], 4:[0,3,9], 5:[0,7,8], 6:[1,8,9], 7:[2,5,9], 8:[3,5,6], 9:[4,6,7]})
            sage.: petersen_spring.show()
            sage: petersen_database = graphs.PetersenGraph()
            sage.: petersen_database.show()
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

################################################################################
#   Families of Graphs
################################################################################

    def CompleteGraph(self, n):
        """
        Returns a complete graph on n nodes.

        A Complete Graph is a graph in which all nodes are connected to all
        other nodes.

        This constructor is dependant on vertices numbered 0 through n-1 in
        NetworkX complete_graph()

        PLOTTING:
        Upon construction, the position dictionary is filled to override
        the spring-layout algorithm. By convention, each complete graph
        will be displayed with the first (0) node at the top, with the
        rest following in a counterclockwise manner.

        In the complete graph, there is a big difference visually in using
        the spring-layout algorithm vs. the position dictionary used in
        this constructor.  The position dictionary flattens the graph,
        making it clear which nodes an edge is connected to.  But the
        complete graph offers a good example of how the spring-layout
        works.  The edges push outward (everything is connected), causing
        the graph to appear as a 3-dimensional pointy ball.  (See examples
        below).

        Filling the position dictionary in advance adds O(n) to the
        constructor.  Feel free to race the constructors below in the
        examples section.  The much larger difference is the time added
        by the spring-layout algorithm when plotting.  (Also shown in the
        example below).  The spring model is typically described as O(n^3),
        as appears to be the case in the NetworkX source code.

        EXAMPLES:
        We view many Complete graphs with a SAGE Graphics Array, first
        with this constructor (i.e., the position dictionary filled):
            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.CompleteGraph(i+3)
            ...    g.append(k)
            ...
            sage.: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(node_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage.: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()

        We compare to plotting with the spring-layout algorithm:
            sage: import networkx
            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    spr = networkx.complete_graph(i+3)
            ...    k = Graph(spr)
            ...    g.append(k)
            ...
            sage.: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(node_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage.: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()
            # Compare the constructors (results will vary)
            sage.: time n = networkx.complete_graph(1559); spring1559 = Graph(n)
            # CPU time: 6.85 s,  Wall time: 9.71 s
            sage.: time posdict1559 = graphs.CompleteGraph(1559)
            #CPU time: 9.67 s,  Wall time: 11.75 s

        We compare the plotting speeds (results will vary):
            sage: n = networkx.complete_graph(23)
            sage: spring23 = Graph(n)
            sage: posdict23 = graphs.CompleteGraph(23)
            sage.: time spring23.show()
            # CPU time: 3.51 s,  Wall time: 4.29 s
            sage.: time posdict23.show()
            # CPU time: 0.82 s,  Wall time: 0.96 s
        """
        pos_dict = {}
        for i in range(n):
            x = float(cos((pi/2) + ((2*pi)/n)*i))
            y = float(sin((pi/2) + ((2*pi)/n)*i))
            pos_dict[i] = [x,y]
        import networkx
        G = networkx.complete_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Complete graph")

    def CompleteBipartiteGraph(self, n1, n2):
        """
        Returns a Complete Bipartite Graph sized n1+n2, with each of the
        nodes [0,(n1-1)] connected to each of the nodes [n1,(n2-1)] and
        vice versa.

        A Complete Bipartite Graph is a graph with its vertices partitioned
        into two groups, V1 and V2.  Each v in V1 is connected to every v
        in V2, and vice versa.

        PLOTTING:
        Upon construction, the position dictionary is filled to override
        the spring-layout algorithm. By convention, each complete bipartite
        graph will be displayed with the first n1 nodes on the top row (at
        y=1) from left to right.  The remaining n2 nodes appear at y=0,
        also from left to right.  The shorter row (partition with fewer
        nodes) is stretched to the same length as the longer row, unless
        the shorter row has 1 node; in which case it is centered.  The x
        values in the plot are in domain [0,max{n1,n2}].

        In the Complete Bipartite graph, there is a visual difference in
        using the spring-layout algorithm vs. the position dictionary used
        in this constructor.  The position dictionary flattens the graph
        and separates the partitioned nodes, making it clear which nodes
        an edge is connected to.  The Complete Bipartite graph plotted with
        the spring-layout algorithm tends to center the nodes in n1 (see
        spring_med in examples below), thus overlapping its nodes and edges,
        making it typically hard to decipher.

        Filling the position dictionary in advance adds O(n) to the
        constructor.  Feel free to race the constructors below in the
        examples section.  The much larger difference is the time added by
        the spring-layout algorithm when plotting.  (Also shown in the
        example below).  The spring model is typically described as O(n^3),
        as appears to be the case in the NetworkX source code.

        EXAMPLES:
        Two ways of constructing the complete bipartite graph, using different
        layout algorithms:
            sage: import networkx
            sage: n = networkx.complete_bipartite_graph(389,157); spring_big = Graph(n)
            sage.: posdict_big = graphs.CompleteBipartiteGraph(389,157)

        Compare the plotting:
            sage: n = networkx.complete_bipartite_graph(11,17)
            sage: spring_med = Graph(n)
            sage: posdict_med = graphs.CompleteBipartiteGraph(11,17)

        Notice here how the spring-layout tends to center the nodes of n1
            sage.: spring_med.show()
            sage.: posdict_med.show()

        View many complete bipartite graphs with a SAGE Graphics Array,
        with this constructor (i.e., the position dictionary filled):
            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.CompleteBipartiteGraph(i+1,4)
            ...    g.append(k)
            ...
            sage.: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(node_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage.: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()

        We compare to plotting with the spring-layout algorithm:
            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    spr = networkx.complete_bipartite_graph(i+1,4)
            ...    k = Graph(spr)
            ...    g.append(k)
            ...
            sage.: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(node_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage.: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()
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
        G = networkx.complete_bipartite_graph(n1,n2)
        return graph.Graph(G, pos=pos_dict, name="Complete bipartite graph")

    def CubeGraph(self, n):
        """
        AUTHOR:  Robert Miller

        PLOTTING:
        See commented source code.

        EXAMPLES:
            # Plot several n-cubes in a SAGE Graphics Array
            sage: g = []
            sage: j = []
            sage: for i in range(6):
            ...    k = graphs.CubeGraph(i+1)
            ...    g.append(k)
            ...
            sage: for i in range(2):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(node_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show(figsize=[6,4])

            # Use the plot options to display larger n-cubes
            sage: g = graphs.CubeGraph(9)
            sage.: g.show(figsize=[12,12],vertex_labels=False, node_size=20)
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

        return graph.Graph(data=d, pos=pos, name="%d-Cube"%n)

    def RandomGNP(self, n, p, seed=None):
        r"""
        Returns a Random graph on $n$ nodes.  Each edge is inserted
        independently with probability $p$.

        IMPLEMENTATION:
        This function calls the NetworkX function \code{gnp_random_graph}.

        REFERENCES:
            P. Erdos and A. Renyi, On Random Graphs, Publ. Math. 6, 290 (1959).
            E. N. Gilbert, Random Graphs, Ann. Math. Stat., 30, 1141 (1959).

        PLOTTING:
        When plotting, this graph will use the default spring-layout
        algorithm, unless a position dictionary is specified.

        EXAMPLES:
        We plot a random graph on 12 nodes with probability $p = .71$:
            sage: gnp = graphs.RandomGNP(12,.71)
            sage.: gnp.show()

        We view many random graphs using a graphics array:
            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.RandomGNP(i+3,.43)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(node_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.plot.GraphicsArray(j)
            sage: G.save('sage.png')

        TIMINGS:
        The following timings compare the speed of RandomGNP and
        RandomGNPFast for sparse and dense graphs:

            time regular_sparse = graphs.RandomGNP(1559,.22)
            # CPU time: 31.79 s,  Wall time: 38.78 s
            time fast_sparse =  graphs.RandomGNPFast(1559,.22)
            # CPU time: 21.72 s,  Wall time: 26.44 s
            time regular_dense = graphs.RandomGNP(1559,.88)
            # CPU time: 38.75 s,  Wall time: 47.65 s
            time fast_dense = graphs.RandomGNP(1559,.88)
            # CPU time: 39.15 s,  Wall time: 48.22 s
        """
        import networkx
        G = networkx.gnp_random_graph(n, p, seed)
        return graph.Graph(G)

    def RandomGNPFast(self, n, p, seed=None):
        """
        Returns a Random graph on $n$ nodes, with each edge inserted
        independently with probability $p$.

        This function calls the NetworkX function \code{fast_gnp_random_graph}.

        PLOTTING:
        When plotting, this graph will use the default spring-layout
        algorithm, unless a position dictionary is specified.

        EXAMPLES:
            # Plot a random graph on 12 nodes with p = .71
            sage: fast = graphs.RandomGNPFast(12,.71)
            sage.: fast.show()

            # View many random graphs using a SAGE Graphics Array
            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.RandomGNPFast(i+3,.43)
            ...    g.append(k)
            ...
            sage.: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(node_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage.: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()

        TIMINGS: See the documentation for RandomGNP.
        """
        import networkx
        G = networkx.fast_gnp_random_graph(n, p, seed)
        return graph.Graph(G)

# Easy access to the graph database from the command line:
graphs = GraphGenerators()





