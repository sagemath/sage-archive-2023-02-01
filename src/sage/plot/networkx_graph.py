from sage.plot.primitive import GraphicPrimitive
from sage.misc.randstate import current_randstate

class NetworkXGraph(GraphicPrimitive):
    """
    Primitive class that initializes the NetworkX graph type.

    INPUT:
        graph -- a NetworkX graph
        pos -- an optional positioning dictionary: for example, the
        spring layout from NetworkX for the 5-cycle is
            {   0: [-0.91679746, 0.88169588,],
                1: [ 0.47294849, 1.125     ,],
                2: [ 1.125     ,-0.12867615,],
                3: [ 0.12743933,-1.125     ,],
                4: [-1.125     ,-0.50118505,]   }
        vertex_labels -- determines whether labels for nodes are plotted
        vertex_size -- node size
        vertex_colors -- a dictionary specifying node colors: each key is a color recognized by
                        matplotlib, and each entry is a list of vertices.
        edge_colors -- a dictionary specifying edge colors: each key is a color recognized by
                        matplotlib, and each entry is a list of edges.
        scaling_term -- default is 0.05. if nodes are getting chopped off, increase; if graph
                        is too small, decrease. should be positive, but values much bigger than
                        1/8 won't be useful unless the nodes are huge
        draw_edges -- whether to draw edges.

    EXAMPLES:
        sage: from sage.plot.networkx_graph import NetworkXGraph
        sage: import networkx
        sage: D = networkx.dodecahedral_graph()
        sage: NGP = NetworkXGraph(D)
        sage: g = Graphics()
        sage: g.add_primitive(NGP)
        sage: g.axes(False)
        sage: g.show()

        sage: import networkx
        sage: from sage.plot.networkx_graph import NetworkXGraph
        sage: from math import sin, cos, pi
        sage: P = networkx.petersen_graph()
        sage: d = {'#FF0000':[0,5], '#FF9900':[1,6], '#FFFF00':[2,7], '#00FF00':[3,8], '#0000FF':[4,9]}
        sage: pos_dict = {}
        sage: for i in range(5):
        ...    x = float(cos(pi/2 + ((2*pi)/5)*i))
        ...    y = float(sin(pi/2 + ((2*pi)/5)*i))
        ...    pos_dict[i] = [x,y]
        ...
        sage: for i in range(10)[5:]:
        ...    x = float(0.5*cos(pi/2 + ((2*pi)/5)*i))
        ...    y = float(0.5*sin(pi/2 + ((2*pi)/5)*i))
        ...    pos_dict[i] = [x,y]
        ...
        sage: NGP = NetworkXGraph(graph=P, vertex_colors=d, pos=pos_dict)
        sage: g = Graphics()
        sage: g.add_primitive(NGP)
        sage: g.axes(False)
        sage: g.show()

        sage: from sage.plot.all import rainbow
        sage: from sage.plot.networkx_graph import NetworkXGraph
        sage: import networkx
        sage: C = graphs.CubeGraph(5)
        sage: pos = C.get_pos()
        sage: G = C.networkx_graph()
        sage: R = rainbow(5)
        sage: edge_colors = {}
        sage: for i in range(5):
        ...    edge_colors[R[i]] = []
        sage: for u,v,l in C.edges():
        ...    for i in range(5):
        ...        if u[i] != v[i]:
        ...            edge_colors[R[i]].append((u,v,l))
        sage: NGP = NetworkXGraph(G, pos=pos, vertex_labels=False, vertex_size=0, edge_colors=edge_colors)
        sage: G = Graphics()
        sage: G.add_primitive(NGP)
        sage: G.axes_range(xmin=-1.1, xmax=2.2, ymin=0, ymax=3.25)
        sage: G.axes(False)
        sage: G.show()

    We color the edges and vertices of a Dodecahedral graph:
        sage: g = graphs.DodecahedralGraph()
        sage: g.show(edge_colors={(1.0, 0.8571428571428571, 0.0): g.edges()})

    """
    def __init__(self, graph, pos=None, vertex_labels=True, vertex_size=300,
                   vertex_colors=None, edge_colors=None, scaling_term=0.05,
                   draw_edges=True):
        self.__nxg = graph
        self.__vertex_size = vertex_size
        self.__vertex_labels = vertex_labels
        self.__vertex_colors = vertex_colors
        self.__edge_colors = edge_colors
        self.__draw_edges = draw_edges
        if len(self.__nxg) != 0:
            import networkx as NX
            if pos is None:
                self.__pos = NX.drawing.spring_layout(self.__nxg)
            else:
                self.__pos = pos

            nodelist=self.__nxg.nodes()

            xes = [self.__pos[v][0] for v in self.__pos]
            ys = [self.__pos[v][1] for v in self.__pos]
            xmin = min(xes)
            xmax = max(xes)
            ymin = min(ys)
            ymax = max(ys)
            if xmax == xmin:
                xmax += 1
                xmin -= 1
            if ymax == ymin:
                ymax += 1
                ymin -= 1
            dx = xmax - xmin
            dy = ymax - ymin

            if not pos is None:
                random = current_randstate().python_random().random
                missing = []
                for v in nodelist:
                    if not v in self.__pos:
                        missing.append(v)
                for v in missing:
                    self.__pos[v] = [xmin + dx*random(),ymin + dy*random()]

            # adjust the plot
            xmin -= scaling_term*dx
            xmax += scaling_term*dx
            ymin -= scaling_term*dy
            ymax += scaling_term*dy

            self._xmin = xmin
            self._xmax = xmax
            self._ymin = ymin
            self._ymax = ymax
        else:
            self.__pos = {}
            self._xmin = -1
            self._xmax = 1
            self._ymin = -1
            self._ymax = 1

    def get_minmax_data(self):
        return {'xmin': self._xmin,
                'xmax': self._xmax,
                'ymin': self._ymin,
                'ymax': self._ymax}

    def _render_on_subplot(self, subplot):
        if len(self.__nxg) != 0:
            import networkx as NX
            vertex_size = float(self.__vertex_size)
            if self.__vertex_colors is None:
                NX.draw_networkx_nodes(G=self.__nxg, pos=self.__pos, ax=subplot, node_size=vertex_size)
            else:
                for i in self.__vertex_colors:
                    NX.draw_networkx_nodes(G=self.__nxg, nodelist=self.__vertex_colors[i],
                                           node_color=i if isinstance(i, str) else [float(z) for z in i],
                                           pos=self.__pos, ax=subplot, node_size=vertex_size)
            if self.__draw_edges:
                if self.__edge_colors is None:
                    NX.draw_networkx_edges(G=self.__nxg, pos=self.__pos, ax=subplot, node_size=vertex_size)
                else:
                    for i in self.__edge_colors:
                        NX.draw_networkx_edges(G=self.__nxg, pos=self.__pos, edgelist=self.__edge_colors[i],
                                               edge_color=i if isinstance(i, str) else [float(z) for z in i],
                                               ax=subplot, node_size=vertex_size)
            if self.__vertex_labels:
                labels = {}
                for v in self.__nxg:
                    labels[v] = str(v)
                NX.draw_networkx_labels(self.__nxg, self.__pos, labels=labels, ax=subplot)

def networkx_plot(graph, pos=None, vertex_labels=True, vertex_size=300, vertex_colors=None,
                  edge_colors=None, graph_border=False, scaling_term=0.05, draw_edges=True):
    """
    Creates a graphics object ready to display a NetworkX graph.

    INPUT:
        graph -- a NetworkX graph
        pos -- an optional positioning dictionary: for example, the
        spring layout from NetworkX for the 5-cycle is
            {   0: [-0.91679746, 0.88169588,],
                1: [ 0.47294849, 1.125     ,],
                2: [ 1.125     ,-0.12867615,],
                3: [ 0.12743933,-1.125     ,],
                4: [-1.125     ,-0.50118505,]   }
        vertex_labels -- determines whether labels for nodes are plotted
        vertex_size -- node size
        vertex_colors -- a dictionary specifying node colors: each key is a color recognized by
                        matplotlib, and each entry is a list of vertices.
        edge_colors -- a dictionary specifying edge colors: each key is a color recognized by
                        matplotlib, and each entry is a list of edges.
        scaling_term -- default is 0.05. if nodes are getting chopped off, increase; if graph
                        is too small, decrease. should be positive, but values much bigger than
                        1/8 won't be useful unless the nodes are huge

    EXAMPLES:
        sage: import networkx
        sage: D = networkx.dodecahedral_graph()
        sage: networkx_plot(D)

        sage: import networkx
        sage: from math import sin, cos, pi
        sage: P = networkx.petersen_graph()
        sage: d = {'#FF0000':[0,5], '#FF9900':[1,6], '#FFFF00':[2,7], '#00FF00':[3,8], '#0000FF':[4,9]}
        sage: pos_dict = {}
        sage: for i in range(5):
        ...    x = float(cos(pi/2 + ((2*pi)/5)*i))
        ...    y = float(sin(pi/2 + ((2*pi)/5)*i))
        ...    pos_dict[i] = [x,y]
        ...
        sage: for i in range(10)[5:]:
        ...    x = float(0.5*cos(pi/2 + ((2*pi)/5)*i))
        ...    y = float(0.5*sin(pi/2 + ((2*pi)/5)*i))
        ...    pos_dict[i] = [x,y]
        ...
        sage: networkx_plot(graph=P, vertex_colors=d, pos=pos_dict)

        sage: C = graphs.CubeGraph(5)
        sage: from sage.plot.plot import rainbow
        sage: R = rainbow(5)
        sage: edge_colors = {}
        sage: for i in range(5):
        ...    edge_colors[R[i]] = []
        sage: for u,v,l in C.edges():
        ...    for i in range(5):
        ...        if u[i] != v[i]:
        ...            edge_colors[R[i]].append((u,v,l))
        sage: networkx_plot(C.networkx_graph(), pos=C.get_pos(), edge_colors=edge_colors, vertex_labels=False, vertex_size=0)

    """
    from sage.plot.plot import Graphics
    g = Graphics()
    NGP = NetworkXGraph(graph, pos=pos, vertex_labels=vertex_labels,
                        vertex_size=vertex_size, vertex_colors=vertex_colors, edge_colors=edge_colors,
                        scaling_term=scaling_term, draw_edges=draw_edges)
    g.add_primitive(NGP)
    xmin = NGP._xmin
    xmax = NGP._xmax
    ymin = NGP._ymin
    ymax = NGP._ymax
    g.axes_range(**NGP.get_minmax_data())
    if graph_border:
        from sage.plot.all import line
        dx = (xmax - xmin)/10
        dy = (ymax - ymin)/10
        border = (line([( xmin - dx, ymin - dy), ( xmin - dx, ymax + dy ), ( xmax + dx, ymax + dy ), ( xmax + dx, ymin - dy ), ( xmin - dx, ymin - dy )], thickness=1.3))
        border.axes_range(xmin = (xmin - dx), xmax = (xmax + dx), ymin = (ymin - dy), ymax = (ymax + dy))
        g = g + border
    g.axes(False)
    return g
