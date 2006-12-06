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

    If you creat a graph in SAGE using the \code{Graph} command, then
    plot that graph, the positioning of nodes is determined using the
    spring-layout algorithm.  For the special graph constructors,
    which you get using \code{graphs.[tab]}, the positions are preset.
    For example, consider the Petersen graph with default node
    positioning vs. the Petersen graph constructed by this database:

        sage: petersen_spring = Graph({0:[1,4,5], 1:[0,2,6], 2:[1,3,7], 3:[2,4,8], 4:[0,3,9], 5:[0,7,8], 6:[1,8,9], 7:[2,5,9], 8:[3,5,6], 9:[4,6,7]})
        sage.: petersen_spring.show()
        sage: petersen_database = graphs.PetersenGraph()
        sage.: petersen_database.show()

    For all the constructors in this database (except the random and
    empty graphs), the position dictionary is filled in, instead of
    using the spring-layout algorithm.

ORGANIZATION:
    The constructors available in this database are organized as follows:
    \begin{verbatim}
        Basic Structures:
            - EmptyGraph
            - CycleGraph
            - StarGraph
            - WheelGraph
        Named Graphs:
            - PetersenGraph
        Families of Graphs:
            - CompleteGraph
            - CompleteBipartiteGraph
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
"""

################################################################################
#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman <eakirkman@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
################################################################################

import networkx
import graph
import sage.functions.functions as functions # sin() and cos()
from   sage.functions.constants import pi

class GraphDatabase():
    r"""
    A class consisting of constructors for several common graphs.

    A list of all graphs and graph structures in this database is available
    via tab completion. Type "graphs." and then hit tab to see which graphs
    are available.

    The docstrings include educational information about each named graph
    with the hopes that this database can be used as a reference.

    For all the constructors in this database (except the random and empty
    graphs), the position dictionary is filled to override the
    spring-layout algorithm.

    The constructors currently in this class include:
    \begin{verbatim}
        Basic Structures:
            - EmptyGraph
            - CycleGraph
            - StarGraph
            - WheelGraph
        Named Graphs:
            - PetersenGraph
        Families of Graphs:
            - CompleteGraph
            - CompleteBipartiteGraph
            - RandomGNP
            - RandomGNPFast
    \end{verbatim}
    """

################################################################################
#   Basic Structures
################################################################################

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
            ...        n.append(g[3*i + m].plot(node_size=50, with_labels=False))
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
            ...        n.append(g[3*i + m].plot(node_size=50, with_labels=False))
            ...    j.append(n)
            ...
            sage.: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()
        """
        pos_dict = {}
        for i in range(n):
            x = float(functions.cos((pi/2) + ((2*pi)/n)*i))
            y = float(functions.sin((pi/2) + ((2*pi)/n)*i))
            pos_dict[i] = [x,y]
        G = networkx.cycle_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Cycle graph on %d vertices"%n)

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
            # Compare the constructors (results will vary)
            sage.: time n = networkx.star_graph(3989); spring3989 = Graph(n)
            # CPU time: 0.08 s,  Wall time: 0.10 s
            sage.: time posdict3989 = graphs.StarGraph(3989)
            # CPU time: 5.43 s,  Wall time: 7.41 s

            # Compare the plotting speeds (results will vary)
            sage: n = networkx.star_graph(23)
            sage: spring23 = Graph(n)
            sage: posdict23 = graphs.StarGraph(23)
            sage.: time spring23.show()
            # CPU time: 2.31 s,  Wall time: 3.14 s
            sage.: time posdict23.show()
            # CPU time: 0.68 s,  Wall time: 0.80 s

            # View many star graphs as a SAGE Graphics Array

            # With this constructor (i.e., the position dictionary filled)
            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.StarGraph(i+3)
            ...    g.append(k)
            ...
            sage.: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(node_size=50, with_labels=False))
            ...    j.append(n)
            ...
            sage.: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()

            # Compared to plotting with the spring-layout algorithm
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
            ...        n.append(g[3*i + m].plot(node_size=50, with_labels=False))
            ...    j.append(n)
            ...
            sage.: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()
        """
        pos_dict = {}
        pos_dict[0] = [0,0]
        for i in range(n+1)[1:]:
            x = float(functions.cos((pi/2) + ((2*pi)/n)*(i-1)))
            y = float(functions.sin((pi/2) + ((2*pi)/n)*(i-1)))
            pos_dict[i] = [x,y]
        G = networkx.star_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Star graph on %d vertices"%(n+1))

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
            ...        n.append(g[3*i + m].plot(node_size=50, with_labels=False))
            ...    j.append(n)
            ...
            sage.: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()

        Next, using the spring-layout algorithm:
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
            ...        n.append(g[3*i + m].plot(node_size=50, with_labels=False))
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
            x = float(functions.cos((pi/2) + ((2*pi)/(n-1))*(i-1)))
            y = float(functions.sin((pi/2) + ((2*pi)/(n-1))*(i-1)))
            pos_dict[i] = [x,y]
        G = networkx.wheel_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Wheel graph on %d vertices"%n)

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
            x = float(functions.cos(pi/2 + ((2*pi)/5)*i))
            y = float(functions.sin(pi/2 + ((2*pi)/5)*i))
            pos_dict[i] = [x,y]
        for i in range(10)[5:]:
            x = float(0.5*functions.cos(pi/2 + ((2*pi)/5)*i))
            y = float(0.5*functions.sin(pi/2 + ((2*pi)/5)*i))
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
            ...        n.append(g[3*i + m].plot(node_size=50, with_labels=False))
            ...    j.append(n)
            ...
            sage.: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()

        We compare to plotting with the spring-layout algorithm:
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
            ...        n.append(g[3*i + m].plot(node_size=50, with_labels=False))
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
            x = float(functions.cos((pi/2) + ((2*pi)/n)*i))
            y = float(functions.sin((pi/2) + ((2*pi)/n)*i))
            pos_dict[i] = [x,y]
        G = networkx.complete_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Complete graph on %d vertices"%n)

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
            ...        n.append(g[3*i + m].plot(node_size=50, with_labels=False))
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
            ...        n.append(g[3*i + m].plot(node_size=50, with_labels=False))
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
        G = networkx.complete_bipartite_graph(n1,n2)
        return graph.Graph(G, pos=pos_dict, name="Complete bipartite graph on %d vertices"%(n1+n2))

    def RandomGNP(self, n, p, seed=None):
        r"""
        Returns a Random graph on $n$ nodes.  Each edge is inserted
        independently with probability $p$.

        IMPLEMENTATION:
        This function calls the NetworkX function \code{gnp_random_graph}.

        REFERENCES:
        \begin{itemize}
           \item P. Erdos and A. Renyi, On Random Graphs, Publ. Math. 6, 290 (1959).
           \item E. N. Gilbert, Random Graphs, Ann. Math. Stat., 30, 1141 (1959).
        \end{itemize}

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
            ...        n.append(g[3*i + m].plot(node_size=50, with_labels=False))
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
            ...        n.append(g[3*i + m].plot(node_size=50, with_labels=False))
            ...    j.append(n)
            ...
            sage.: G = sage.plot.plot.GraphicsArray(j)
            sage.: G.show()

        TIMINGS: See the documentation for RandomGNP.
        """
        G = networkx.fast_gnp_random_graph(n, p, seed)
        return graph.Graph(G)

# Easy access to the graph database from the command line:
graphs = GraphDatabase()





