"""
Graph Plotting
"""

#*****************************************************************************
#      Copyright (C) 2009   Emily Kirkman
#                    2009   Robert L. Miller <rlmillster@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.sage_object import SageObject
from sage.plot.all import Graphics, scatter_plot, bezier_path, line, arrow, text, circle
from sage.misc.decorators import options
from math import sqrt, cos, sin, atan, pi

layout_options =   {
                    'layout': 'A layout algorithm -- one of "acyclic", "circular", "ranked", "graphviz", "planar", "spring", or "tree".',
                    'iterations': 'The number of times to execute the spring layout algorithm.',
                    'heights': 'A dictionary mapping heights to the list of vertices at this height.',
                    'spring': 'Use spring layout to finalize the current layout.',
                    'tree_root': 'A vertex designation for drawing trees.',
                    'tree_orientation': 'The direction of tree branches -- "up" or "down".',
                    'save_pos': 'Whether or not to save the computed position for the graph.',
                    'dim': 'The dimension of the layout -- 2 or 3.',
                    'prog': 'Which graphviz layout program to use -- one of "circo", "dot", "fdp", "neato", or "twopi".',
                    'by_component': 'Whether to do the spring layout by connected component -- a boolean.',
                    }

graphplot_options = layout_options.copy()
graphplot_options.update(
                   {'pos': 'The position dictionary of vertices',
                    'vertex_labels': 'Whether or not to draw vertex labels.',
                    'vertex_colors': 'Dictionary of vertex coloring.',
                    'vertex_size': 'The size to draw the vertices.',
                    'vertex_shape': 'The shape to draw the vertices, Currently unavailable for Multi-edged DiGraphs.',
                    'edge_labels': 'Whether or not to draw edge labels.',
                    'edge_style': 'The linestyle of the edges-- one of "solid", "dashed", "dotted", dashdot".',
                    'edge_color': 'The default color for edges.',
                    'edge_colors': 'Dictionary of edge coloring.',
                    'color_by_label': 'Whether or not to color the edges by their label values.',
                    'partition': 'A partition of the vertex set.  (Draws each cell of vertices in a different color).',
                    'loop_size': 'The radius of the smallest loop.',
                    'dist': 'The distance between multiedges.',
                    'max_dist': 'The max distance range to allow multiedges.',
                    'talk': 'Whether to display the vertices in talk mode (larger and white)',
                    'graph_border': 'Whether or not to draw a frame around the graph.'})

class GraphPlot(SageObject):
    def __init__(self, graph, options):
        """
        Returns a GraphPlot object, which stores all the parameters needed for
        plotting (Di)Graphs.  A GraphPlot has a plot and show function, as well
        as some functions to set parameters for vertices and edges.  This constructor
        assumes default options are set.  Defaults are shown in the example below.

        EXAMPLE::

            sage: from sage.graphs.graph_plot import GraphPlot
            sage: options = {
            ...     'vertex_size':200,
            ...     'vertex_labels':True,
            ...     'layout':None,
            ...     'edge_style':'solid',
            ...     'edge_color':'black',
            ...     'edge_colors':None,
            ...     'edge_labels':False,
            ...     'iterations':50,
            ...     'tree_orientation':'down',
            ...     'heights':None,
            ...     'graph_border':False,
            ...     'talk':False,
            ...     'color_by_label':False,
            ...     'partition':None,
            ...     'dist':.075,
            ...     'max_dist':1.5,
            ...     'loop_size':.075}
            sage: g = Graph({0:[1,2], 2:[3], 4:[0,1]})
            sage: GP = GraphPlot(g, options)

        TESTS::

            sage: g = graphs.CompleteGraph(2); g.show()

        """
        self._plot_components = {}
        self._nodelist = graph.vertices()
        self._graph = graph
        self._options = options
        self.set_pos()
        self._arcs = self._graph.has_multiple_edges(to_undirected=True)
        self._loops = self._graph.has_loops()
        if self._graph.is_directed() and self._arcs:
            self._arcdigraph = True
        else:
            self._arcdigraph = False
        self.set_vertices()
        self.set_edges()

    def _repr_(self):
        """
        Returns a string representation of a GraphPlot object.

        EXAMPLE:

        This function is called implicitly by the code below::

            sage: g = Graph({0:[1,2], 2:[3], 4:[0,1]})
            sage: g.graphplot()
            GraphPlot object for Graph on 5 vertices
        """
        return "GraphPlot object for %s"%self._graph

    def set_pos(self):
        """
        Sets the position plotting parameters for this GraphPlot.

        EXAMPLES:

        This function is called implicitly by the code below::

            sage: g = Graph({0:[1,2], 2:[3], 4:[0,1]})
            sage: g.graphplot(save_pos=True, layout='circular')
            GraphPlot object for Graph on 5 vertices

        The following illustrates the format of a position dictionary,
        but due to numerical noise we do not check the values themselves.

        ::

            sage: g.get_pos()
            {0: [...e-17, 1.0],
             1: [-0.951..., 0.309...],
             2: [-0.587..., -0.809...],
             3: [0.587..., -0.809...],
             4: [0.951..., 0.309...]}

            sage: T = list(graphs.trees(7))
            sage: t = T[3]
            sage: t.plot(heights={0:[0], 1:[4,5,1], 2:[2], 3:[3,6]})

        TESTS:

        Make sure that vertex locations are floats.  Not being floats
        isn't a bug in itself but makes it too easy to accidentally
        introduce a bug elsewhere, such as in set_edges (trac #10124),
        via silent truncating division of integers::

            sage: g = graphs.FruchtGraph()
            sage: gp = g.graphplot()
            sage: set(map(type, flatten(gp._pos.values())))
            set([<type 'float'>])
            sage: g = graphs.BullGraph()
            sage: gp = g.graphplot(save_pos=True)
            sage: set(map(type, flatten(gp._pos.values())))
            set([<type 'float'>])

        """
        self._pos = self._graph.layout(**self._options)
        # make sure the positions are floats (trac #10124)
        self._pos = dict((k,(float(v[0]), float(v[1]))) for k,v in self._pos.iteritems())

    def set_vertices(self, **vertex_options):
        """
        Sets the vertex plotting parameters for this GraphPlot.  This function
        is called by the constructor but can also be called to make updates to
        the vertex options of an existing GraphPlot object.  Note that the
        changes are cumulative.

        EXAMPLES::

            sage: g = Graph({}, loops=True, multiedges=True, sparse=True)
            sage: g.add_edges([(0,0,'a'),(0,0,'b'),(0,1,'c'),(0,1,'d'),
            ...     (0,1,'e'),(0,1,'f'),(0,1,'f'),(2,1,'g'),(2,2,'h')])
            sage: GP = g.graphplot(vertex_size=100, edge_labels=True, color_by_label=True, edge_style='dashed')
            sage: GP.set_vertices(talk=True)
            sage: GP.plot()
            sage: GP.set_vertices(vertex_colors='pink', vertex_shape='^')
            sage: GP.plot()
        """
        # Handle base vertex options
        voptions = {}

        for arg in vertex_options:
            self._options[arg] = vertex_options[arg]

        # First set defaults for styles
        vertex_colors = None
        if self._options['talk']:
            voptions['markersize'] = 500
            if self._options['partition'] is None:
                vertex_colors = '#ffffff'
        else:
            voptions['markersize'] = self._options['vertex_size']

        if 'vertex_colors' not in self._options or self._options['vertex_colors'] is None:
            if self._options['partition'] is not None:
                from sage.plot.colors import rainbow,rgbcolor
                partition = self._options['partition']
                l = len(partition)
                R = rainbow(l)
                vertex_colors = {}
                for i in range(l):
                    vertex_colors[R[i]] = partition[i]
            elif len(self._graph._boundary) != 0:
                vertex_colors = {}
                bdy_verts = []
                int_verts = []
                for v in self._graph.vertex_iterator():
                    if v in self._graph._boundary:
                        bdy_verts.append(v)
                    else:
                        int_verts.append(v)
                vertex_colors['#fec7b8'] = int_verts
                vertex_colors['#b3e8ff'] = bdy_verts
            elif not vertex_colors:
                vertex_colors='#fec7b8'
        else:
            vertex_colors = self._options['vertex_colors']

        if 'vertex_shape' in self._options:
            voptions['marker'] = self._options['vertex_shape']

        if self._graph.is_directed():
            self._vertex_radius = sqrt(voptions['markersize']/pi)
            self._arrowshorten = 2*self._vertex_radius
            if self._arcdigraph:
                self._vertex_radius = sqrt(voptions['markersize']/(20500*pi))

        voptions['zorder'] = 7

        if not isinstance(vertex_colors, dict):
            voptions['facecolor'] = vertex_colors
            if self._arcdigraph:
                self._plot_components['vertices'] = [circle(center,
                    self._vertex_radius, fill=True, facecolor=vertex_colors, clip=False)
                    for center in self._pos.values()]
            else:
                self._plot_components['vertices'] = scatter_plot(
                    self._pos.values(), clip=False, **voptions)
        else:
            # Color list must be ordered:
            pos = []
            colors = []
            for i in vertex_colors:
                pos += [self._pos[j] for j in vertex_colors[i]]
                colors += [i]*len(vertex_colors[i])

            # If all the vertices have not been assigned a color
            if len(self._pos)!=len(pos):
                from sage.plot.colors import rainbow,rgbcolor
                vertex_colors_rgb=[rgbcolor(c) for c in vertex_colors]
                for c in rainbow(len(vertex_colors)+1):
                    if rgbcolor(c) not in vertex_colors_rgb:
                        break
                leftovers=[j for j in self._pos.values() if j not in pos]
                pos+=leftovers
                colors+=[c]*len(leftovers)

            if self._arcdigraph:
                self._plot_components['vertices'] = [circle(pos[i],
                    self._vertex_radius, fill=True, facecolor=colors[i], clip=False)
                    for i in range(len(pos))]
            else:
                self._plot_components['vertices'] = scatter_plot(pos,
                    facecolor=colors, clip=False, **voptions)

        if self._options['vertex_labels']:
            self._plot_components['vertex_labels'] = []
            # TODO: allow text options
            for v in self._nodelist:
                self._plot_components['vertex_labels'].append(text(str(v),
                    self._pos[v], rgbcolor=(0,0,0), zorder=8))

    def set_edges(self, **edge_options):
        """
        Sets the edge (or arrow) plotting parameters for the GraphPlot object.  This
        function is called by the constructor but can also be called to make updates to
        the vertex options of an existing GraphPlot object.  Note that the changes are
        cumulative.

        EXAMPLES::

            sage: g = Graph({}, loops=True, multiedges=True, sparse=True)
            sage: g.add_edges([(0,0,'a'),(0,0,'b'),(0,1,'c'),(0,1,'d'),
            ...     (0,1,'e'),(0,1,'f'),(0,1,'f'),(2,1,'g'),(2,2,'h')])
            sage: GP = g.graphplot(vertex_size=100, edge_labels=True, color_by_label=True, edge_style='dashed')
            sage: GP.set_edges(edge_style='solid')
            sage: GP.plot()
            sage: GP.set_edges(edge_color='black')
            sage: GP.plot()

            sage: d = DiGraph({}, loops=True, multiedges=True, sparse=True)
            sage: d.add_edges([(0,0,'a'),(0,0,'b'),(0,1,'c'),(0,1,'d'),
            ...     (0,1,'e'),(0,1,'f'),(0,1,'f'),(2,1,'g'),(2,2,'h')])
            sage: GP = d.graphplot(vertex_size=100, edge_labels=True, color_by_label=True, edge_style='dashed')
            sage: GP.set_edges(edge_style='solid')
            sage: GP.plot()
            sage: GP.set_edges(edge_color='black')
            sage: GP.plot()

        TESTS::

            sage: G = Graph("Fooba")
            sage: G.show(edge_colors={'red':[(3,6),(2,5)]})

        Verify that default edge labels are pretty close to being between the vertices
        in some cases where they weren't due to truncating division (trac #10124)::

            sage: test_graphs = graphs.FruchtGraph(), graphs.BullGraph()
            sage: tol = 0.001
            sage: for G in test_graphs:
            ...       E=G.edges()
            ...       for e0, e1, elab in E:
            ...           G.set_edge_label(e0, e1, '%d %d' % (e0, e1))
            ...       gp = G.graphplot(save_pos=True,edge_labels=True)
            ...       vx = gp._plot_components['vertices'][0].xdata
            ...       vy = gp._plot_components['vertices'][0].ydata
            ...       for elab in gp._plot_components['edge_labels']:
            ...           textobj = elab[0]
            ...           x, y, s = textobj.x, textobj.y, textobj.string
            ...           v0, v1 = map(int, s.split())
            ...           vn = vector(((x-(vx[v0]+vx[v1])/2.),y-(vy[v0]+vy[v1])/2.)).norm()
            ...           assert vn < tol


        """
        for arg in edge_options:
            self._options[arg] = edge_options[arg]
        if 'edge_colors' in edge_options: self._options['color_by_label'] = False

        # Handle base edge options: thickness, linestyle
        eoptions={}
        if 'edge_style' in self._options:
            eoptions['linestyle'] = self._options['edge_style']
        if 'thickness' in self._options:
            eoptions['thickness'] = self._options['thickness']

        # Set labels param to add labels on the fly
        labels = False
        if self._options['edge_labels']:
            labels = True
            self._plot_components['edge_labels'] = []

        # Make dict collection of all edges (keep label and edge color)
        edges_to_draw = {}
        if self._options['color_by_label'] or isinstance(self._options['edge_colors'], dict):
            if self._options['color_by_label']: edge_colors = self._graph._color_by_label()
            else: edge_colors = self._options['edge_colors']
            for color in edge_colors:
                for edge in edge_colors[color]:
                    key = tuple(sorted([edge[0],edge[1]]))
                    if key == (edge[0],edge[1]): head = 1
                    else: head = 0

                    if len(edge) < 3:
                        label = self._graph.edge_label(edge[0],edge[1])
                        if isinstance(label, list):
                            if key in edges_to_draw:
                                edges_to_draw[key].append((label[-1], color, head))
                            else:
                                edges_to_draw[key] = [(label[-1], color, head)]
                            for i in range(len(label)-1):
                                edges_to_draw[key].append((label[-1], color, head))
                    else:
                        label = edge[2]

                    if key in edges_to_draw:
                        edges_to_draw[key].append((label, color, head))
                    else:
                        edges_to_draw[key] = [(label, color, head)]
            # add unspecified edges in (default color black)
            for edge in self._graph.edge_iterator():
                key = tuple(sorted([edge[0],edge[1]]))
                label = edge[2]
                specified = False
                if key in edges_to_draw:
                    for old_label, old_color, old_head in edges_to_draw[key]:
                        if label == old_label:
                            specified = True
                            break
                if not specified:
                    if key == (edge[0],edge[1]): head = 1
                    else: head = 0
                    edges_to_draw[key] = [(label, 'black', head)]
        else:
            for edge in self._graph.edges(sort=True):
                key = tuple(sorted([edge[0],edge[1]]))
                if key == (edge[0],edge[1]): head = 1
                else: head = 0
                if key in edges_to_draw:
                    edges_to_draw[key].append((edge[2], self._options['edge_color'], head))
                else:
                    edges_to_draw[key] = [(edge[2], self._options['edge_color'], head)]

        if edges_to_draw:
            self._plot_components['edges'] = []
        else:
            return

        # Check for multi-edges or loops
        if self._arcs or self._loops:
            tmp = edges_to_draw.copy()
            dist = self._options['dist']*2.
            loop_size = self._options['loop_size']
            max_dist = self._options['max_dist']
            from sage.functions.all import sqrt
            for (a,b) in tmp:
                if a == b:
                    # Loops
                    distance = dist
                    local_labels = edges_to_draw.pop((a,b))
                    if len(local_labels)*dist > max_dist:
                        distance = float(max_dist)/len(local_labels)
                    curr_loop_size = loop_size
                    for i in range(len(local_labels)):
                        self._plot_components['edges'].append(circle((self._pos[a][0],self._pos[a][1]-curr_loop_size), curr_loop_size, rgbcolor=local_labels[i][1], **eoptions))
                        if labels:
                            self._plot_components['edge_labels'].append(text(local_labels[i][0], (self._pos[a][0], self._pos[a][1]-2*curr_loop_size)))
                        curr_loop_size += distance/4
                elif len(edges_to_draw[(a,b)]) > 1:
                    # Multi-edge
                    local_labels = edges_to_draw.pop((a,b))

                    # Compute perpendicular bisector
                    p1 = self._pos[a]
                    p2 = self._pos[b]
                    M = ((p1[0]+p2[0])/2., (p1[1]+p2[1])/2.) # midpoint
                    if not p1[1] == p2[1]:
                        S = float(p1[0]-p2[0])/(p2[1]-p1[1]) # perp slope
                        y = lambda x : S*x-S*M[0]+M[1] # perp bisector line

                        # f,g are functions of distance d to determine x values
                        # on line y at d from point M
                        f = lambda d : sqrt(d**2/(1.+S**2)) + M[0]
                        g = lambda d : -sqrt(d**2/(1.+S**2)) + M[0]

                        odd_x = f
                        even_x = g
                        if p1[0] == p2[0]:
                            odd_y = lambda d : M[1]
                            even_y = odd_y
                        else:
                            odd_y = lambda x : y(f(x))
                            even_y = lambda x : y(g(x))
                    else:
                        odd_x = lambda d : M[0]
                        even_x = odd_x
                        odd_y = lambda d : M[1] + d
                        even_y = lambda d : M[1] - d

                    # We now have the control points for each bezier curve
                    # in terms of distance parameter d.
                    # Also note that the label for each edge should be drawn at d/2.
                    # (This is because we're using the perp bisectors).
                    distance = dist
                    if len(local_labels)*dist > max_dist:
                        distance = float(max_dist)/len(local_labels)
                    for i in range(len(local_labels)/2):
                        k = (i+1.0)*distance
                        if self._arcdigraph:
                            odd_start = self._polar_hack_for_multidigraph(p1, [odd_x(k),odd_y(k)], self._vertex_radius)[0]
                            odd_end = self._polar_hack_for_multidigraph([odd_x(k),odd_y(k)], p2, self._vertex_radius)[1]
                            even_start = self._polar_hack_for_multidigraph(p1, [even_x(k),even_y(k)], self._vertex_radius)[0]
                            even_end = self._polar_hack_for_multidigraph([even_x(k),even_y(k)], p2, self._vertex_radius)[1]
                            self._plot_components['edges'].append(arrow(path=[[odd_start,[odd_x(k),odd_y(k)],odd_end]], head=local_labels[2*i][2], zorder=1, rgbcolor=local_labels[2*i][1], **eoptions))
                            self._plot_components['edges'].append(arrow(path=[[even_start,[even_x(k),even_y(k)],even_end]], head=local_labels[2*i+1][2], zorder=1, rgbcolor=local_labels[2*i+1][1], **eoptions))
                        else:
                            self._plot_components['edges'].append(bezier_path([[p1,[odd_x(k),odd_y(k)],p2]],zorder=1, rgbcolor=local_labels[2*i][1], **eoptions))
                            self._plot_components['edges'].append(bezier_path([[p1,[even_x(k),even_y(k)],p2]],zorder=1, rgbcolor=local_labels[2*i+1][1], **eoptions))
                        if labels:
                            j = k/2.0
                            self._plot_components['edge_labels'].append(text(local_labels[2*i][0],[odd_x(j),odd_y(j)]))
                            self._plot_components['edge_labels'].append(text(local_labels[2*i+1][0],[even_x(j),even_y(j)]))
                    if len(local_labels)%2 == 1:
                        edges_to_draw[(a,b)] = [local_labels[-1]] # draw line for last odd

        dir = self._graph.is_directed()
        for (a,b) in edges_to_draw:
            if self._arcdigraph:
                C,D = self._polar_hack_for_multidigraph(self._pos[a], self._pos[b], self._vertex_radius)
                self._plot_components['edges'].append(arrow(C,D, rgbcolor=edges_to_draw[(a,b)][0][1], head=edges_to_draw[(a,b)][0][2], **eoptions))
                if labels:
                    self._plot_components['edge_labels'].append(text(str(edges_to_draw[(a,b)][0][0]),[(C[0]+D[0])/2., (C[1]+D[1])/2.]))
            elif dir:
                self._plot_components['edges'].append(arrow(self._pos[a],self._pos[b], rgbcolor=edges_to_draw[(a,b)][0][1], arrowshorten=self._arrowshorten, head=edges_to_draw[(a,b)][0][2], **eoptions))
            else:
                self._plot_components['edges'].append(line([self._pos[a],self._pos[b]], rgbcolor=edges_to_draw[(a,b)][0][1], **eoptions))
            if labels and not self._arcdigraph:
                self._plot_components['edge_labels'].append(text(str(edges_to_draw[(a,b)][0][0]),[(self._pos[a][0]+self._pos[b][0])/2., (self._pos[a][1]+self._pos[b][1])/2.]))

    def _polar_hack_for_multidigraph(self, A, B, VR):
        """
        Helper function to quickly compute the two points of intersection of a line
        segment from A to B (xy tuples) and circles centered at A and B, both with
        radius VR.  Returns a tuple of xy tuples representing the two points.

        EXAMPLE::

            sage: d = DiGraph({}, loops=True, multiedges=True, sparse=True)
            sage: d.add_edges([(0,0,'a'),(0,0,'b'),(0,1,'c'),(0,1,'d'),
            ...     (0,1,'e'),(0,1,'f'),(0,1,'f'),(2,1,'g'),(2,2,'h')])
            sage: GP = d.graphplot(vertex_size=100, edge_labels=True, color_by_label=True, edge_style='dashed')
            sage: GP._polar_hack_for_multidigraph((0,1),(1,1),.1)
            ([0.10..., 1.00...], [0.90..., 1.00...])

        TESTS:

        Make sure that Python ints are acceptable arguments (trac #10124)::

            sage: GP = DiGraph().graphplot()
            sage: GP._polar_hack_for_multidigraph((0, 1), (2, 2), .1)
            ([0.08..., 1.04...], [1.91..., 1.95...])
            sage: GP._polar_hack_for_multidigraph((int(0),int(1)),(int(2),int(2)),.1)
            ([0.08..., 1.04...], [1.91..., 1.95...])

        """
        D = [float(B[i]-A[i]) for i in range(2)]
        R = sqrt(D[0]**2+D[1]**2)
        theta = 3*pi/2
        if D[0] > 0:
            theta = atan(D[1]/D[0])
            if D[1] < 0:
                theta += 2*pi
        elif D[0] < 0:
            theta = atan(D[1]/D[0]) + pi
        elif D[1] > 0:
            theta = pi/2
        return ([VR*cos(theta)+A[0], VR*sin(theta)+A[1]], [(R-VR)*cos(theta)+A[0], (R-VR)*sin(theta)+A[1]])

    def show(self, **kwds):
        """
        Shows the (Di)Graph associated with this GraphPlot object.

        For syntax and lengthy documentation, see GP.plot?. Any options not used by
        plot will be passed on to the Graphics.show method.

        EXAMPLE::

            sage: C = graphs.CubeGraph(8)
            sage: P = C.graphplot(vertex_labels=False, vertex_size=0, graph_border=True)
            sage: P.show()
        """
        self.plot().show(**kwds)

    def plot(self, **kwds):
        """
        Returns a graphics object representing the (di)graph.

        INPUT:
            - pos -- an optional positioning dictionary
            - layout -- what kind of layout to use, takes precedence over pos

              - 'circular' -- plots the graph with vertices evenly distributed
                on a circle
              - 'spring' -- uses the traditional spring layout, using the
                graph's current positions as initial positions
              - 'tree' -- the (di)graph must be a tree. One can specify the root
                of the tree using the keyword tree_root, otherwise a root
                will be selected at random. Then the tree will be plotted in
                levels, depending on minimum distance for the root.
            - vertex_labels -- whether to print vertex labels
              edge_labels -- whether to print edge labels. By default, False,
              but if True, the result of str(l) is printed on the edge for
              each label l. Labels equal to None are not printed (to set edge
              labels, see set_edge_label).
            - vertex_size -- size of vertices displayed
            - vertex_shape -- the shape to draw the vertices (Not available for
              multiedge digraphs.
            - graph_border -- whether to include a box around the graph
            - vertex_colors -- optional dictionary to specify vertex colors: each
              key is a color recognizable by matplotlib, and each corresponding
              entry is a list of vertices. If a vertex is not listed, it looks
              invisible on the resulting plot (it doesn't get drawn).
            - edge_colors -- a dictionary specifying edge colors: each key is a
              color recognized by matplotlib, and each entry is a list of edges.
            - partition -- a partition of the vertex set. if specified, plot will
              show each cell in a different color. vertex_colors takes precedence.
            - talk -- if true, prints large vertices with white backgrounds so that
              labels are legible on slides
            - iterations -- how many iterations of the spring layout algorithm to
              go through, if applicable
            - color_by_label -- if True, color edges by their labels
            - heights -- if specified, this is a dictionary from a set of
              floating point heights to a set of vertices
            - edge_style -- keyword arguments passed into the
              edge-drawing routine.  This currently only works for
              directed graphs, since we pass off the undirected graph to
              networkx
            - tree_root -- a vertex of the tree to be used as the root for
              the layout="tree" option. If no root is specified, then one
              is chosen at random. Ignored unless layout='tree'.
            - tree_orientation -- "up" or "down" (default is "down").
              If "up" (resp., "down"), then the root of the tree will
              appear on the bottom (resp., top) and the tree will grow
              upwards (resp. downwards). Ignored unless layout='tree'.
            - save_pos -- save position computed during plotting

        EXAMPLES::

            sage: from sage.graphs.graph_plot import graphplot_options
            sage: list(sorted(graphplot_options.iteritems()))
            [('by_component', 'Whether to do the spring layout by connected component -- a boolean.'),
             ('color_by_label', 'Whether or not to color the edges by their label values.'),
             ('dim', 'The dimension of the layout -- 2 or 3.'),
             ('dist', 'The distance between multiedges.'),
             ('edge_color', 'The default color for edges.'),
             ('edge_colors', 'Dictionary of edge coloring.'),
             ('edge_labels', 'Whether or not to draw edge labels.'),
             ('edge_style', 'The linestyle of the edges-- one of "solid", "dashed", "dotted", dashdot".'),
             ('graph_border', 'Whether or not to draw a frame around the graph.'),
             ('heights', 'A dictionary mapping heights to the list of vertices at this height.'),
             ('iterations', 'The number of times to execute the spring layout algorithm.'),
             ('layout', 'A layout algorithm -- one of "acyclic", "circular", "ranked", "graphviz", "planar", "spring", or "tree".'),
             ('loop_size', 'The radius of the smallest loop.'),
             ('max_dist', 'The max distance range to allow multiedges.'),
             ('partition', 'A partition of the vertex set.  (Draws each cell of vertices in a different color).'),
             ('pos', 'The position dictionary of vertices'),
             ('prog', 'Which graphviz layout program to use -- one of "circo", "dot", "fdp", "neato", or "twopi".'),
             ('save_pos', 'Whether or not to save the computed position for the graph.'),
             ('spring', 'Use spring layout to finalize the current layout.'),
             ('talk', 'Whether to display the vertices in talk mode (larger and white)'),
             ('tree_orientation', 'The direction of tree branches -- "up" or "down".'),
             ('tree_root', 'A vertex designation for drawing trees.'),
             ('vertex_colors', 'Dictionary of vertex coloring.'),
             ('vertex_labels', 'Whether or not to draw vertex labels.'),
             ('vertex_shape', 'The shape to draw the vertices, Currently unavailable for Multi-edged DiGraphs.'),
             ('vertex_size', 'The size to draw the vertices.')]

            sage: from math import sin, cos, pi
            sage: P = graphs.PetersenGraph()
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
            sage: pl = P.graphplot(pos=pos_dict, vertex_colors=d)
            sage: pl.show()

            sage: C = graphs.CubeGraph(8)
            sage: P = C.graphplot(vertex_labels=False, vertex_size=0, graph_border=True)
            sage: P.show()

            sage: G = graphs.HeawoodGraph().copy(sparse=True)
            sage: for u,v,l in G.edges():
            ...    G.set_edge_label(u,v,'(' + str(u) + ',' + str(v) + ')')
            sage: G.graphplot(edge_labels=True).show()

            sage: D = DiGraph( { 0: [1, 10, 19], 1: [8, 2], 2: [3, 6], 3: [19, 4], 4: [17, 5], 5: [6, 15], 6: [7], 7: [8, 14], 8: [9], 9: [10, 13], 10: [11], 11: [12, 18], 12: [16, 13], 13: [14], 14: [15], 15: [16], 16: [17], 17: [18], 18: [19], 19: []}, implementation='networkx' )
            sage: for u,v,l in D.edges():
            ...    D.set_edge_label(u,v,'(' + str(u) + ',' + str(v) + ')')
            sage: D.graphplot(edge_labels=True, layout='circular').show()

            sage: from sage.plot.colors import rainbow
            sage: C = graphs.CubeGraph(5)
            sage: R = rainbow(5)
            sage: edge_colors = {}
            sage: for i in range(5):
            ...    edge_colors[R[i]] = []
            sage: for u,v,l in C.edges():
            ...    for i in range(5):
            ...        if u[i] != v[i]:
            ...            edge_colors[R[i]].append((u,v,l))
            sage: C.graphplot(vertex_labels=False, vertex_size=0, edge_colors=edge_colors).show()

            sage: D = graphs.DodecahedralGraph()
            sage: Pi = [[6,5,15,14,7],[16,13,8,2,4],[12,17,9,3,1],[0,19,18,10,11]]
            sage: D.show(partition=Pi)

            sage: G = graphs.PetersenGraph()
            sage: G.allow_loops(True)
            sage: G.add_edge(0,0)
            sage: G.show()

            sage: D = DiGraph({0:[0,1], 1:[2], 2:[3]}, loops=True)
            sage: D.show()
            sage: D.show(edge_colors={(0,1,0):[(0,1,None),(1,2,None)],(0,0,0):[(2,3,None)]})

            sage: pos = {0:[0.0, 1.5], 1:[-0.8, 0.3], 2:[-0.6, -0.8], 3:[0.6, -0.8], 4:[0.8, 0.3]}
            sage: g = Graph({0:[1], 1:[2], 2:[3], 3:[4], 4:[0]})
            sage: g.graphplot(pos=pos, layout='spring', iterations=0).plot()

            sage: G = Graph()
            sage: P = G.graphplot().plot()
            sage: P.axes()
            False
            sage: G = DiGraph()
            sage: P = G.graphplot().plot()
            sage: P.axes()
            False

            sage: T = list(graphs.trees(7))
            sage: t = T[3]
            sage: t.graphplot(heights={0:[0], 1:[4,5,1], 2:[2], 3:[3,6]}).plot()

            sage: T = list(graphs.trees(7))
            sage: t = T[3]
            sage: t.graphplot(heights={0:[0], 1:[4,5,1], 2:[2], 3:[3,6]}).plot()
            sage: t.set_edge_label(0,1,-7)
            sage: t.set_edge_label(0,5,3)
            sage: t.set_edge_label(0,5,99)
            sage: t.set_edge_label(1,2,1000)
            sage: t.set_edge_label(3,2,'spam')
            sage: t.set_edge_label(2,6,3/2)
            sage: t.set_edge_label(0,4,66)
            sage: t.graphplot(heights={0:[0], 1:[4,5,1], 2:[2], 3:[3,6]}, edge_labels=True).plot()

            sage: T = list(graphs.trees(7))
            sage: t = T[3]
            sage: t.graphplot(layout='tree').show()

            sage: t = DiGraph('JCC???@A??GO??CO??GO??')
            sage: t.graphplot(layout='tree', tree_root=0, tree_orientation="up").show()

            sage: D = DiGraph({0:[1,2,3], 2:[1,4], 3:[0]})
            sage: D.graphplot().show()

            sage: D = DiGraph(multiedges=True, sparse=True)
            sage: for i in range(5):
            ...     D.add_edge((i,i+1,'a'))
            ...     D.add_edge((i,i-1,'b'))
            sage: D.graphplot(edge_labels=True,edge_colors=D._color_by_label()).plot()

            sage: g = Graph({}, loops=True, multiedges=True, sparse=True)
            sage: g.add_edges([(0,0,'a'),(0,0,'b'),(0,1,'c'),(0,1,'d'),
            ...     (0,1,'e'),(0,1,'f'),(0,1,'f'),(2,1,'g'),(2,2,'h')])
            sage: g.graphplot(edge_labels=True, color_by_label=True, edge_style='dashed').plot()
        """
        G = Graphics()
        for comp in self._plot_components.values():
            if not isinstance(comp, list):
                G += comp
            else:
                for item in comp:
                    G += item
        G.set_axes_range(*(self._graph._layout_bounding_box(self._pos)))
        if self._options['graph_border']:
            xmin = G.xmin()
            xmax = G.xmax()
            ymin = G.ymin()
            ymax = G.ymax()
            dx = (xmax-xmin)/10.0
            dy = (ymax-ymin)/10.0
            border = (line([( xmin - dx, ymin - dy), ( xmin - dx, ymax + dy ), ( xmax + dx, ymax + dy ), ( xmax + dx, ymin - dy ), ( xmin - dx, ymin - dy )], thickness=1.3))
            border.axes_range(xmin = (xmin - dx), xmax = (xmax + dx), ymin = (ymin - dy), ymax = (ymax + dy))
            G += border
        G.set_aspect_ratio(1)
        G.axes(False)
        G._extra_kwds['axes_pad']=.05
        return G

    def layout_tree(self,root,orientation):
        """

        Compute a nice layout of a tree.

        INPUT:

        - ``root`` -- the root vertex.

        - ``orientation`` -- Whether to place the root
          at the top or at the bottom :

            - ``orientation="down"`` -- children are placed below
              their parent
            - ``orientation="top"`` -- children are placed above
              their parent


        EXAMPLES::

            sage: T = graphs.RandomLobster(25,0.3,0.3)
            sage: T.show(layout='tree',tree_orientation='up')

            sage: from sage.graphs.graph_plot import GraphPlot
            sage: G = graphs.HoffmanSingletonGraph()
            sage: T = Graph()
            sage: T.add_edges(G.min_spanning_tree(starting_vertex=0))
            sage: T.show(layout='tree',tree_root=0)

        """

        T = self._graph

        if not self._graph.is_tree():
            raise RuntimeError("Cannot use tree layout on this graph: self.is_tree() returns False.")

        children = {root:T.neighbors(root)}

        #always make a copy of the children because they get eaten
        stack = [[u for u in children[root]]]
        stick = [root]
        parent = dict([(u,root) for u in children[root]])
        pos = {}
        obstruction = [0.0]*T.num_verts()
        if orientation == 'down':
            o = -1
        else:
            o = 1

        def slide(v,dx):
            """

            shift the vertex v and its descendants to the right by dx

            Precondition: v and its descendents have already had their
            positions computed.

            """

            level = [v]
            while level:
                nextlevel = []
                for u in level:
                    x,y = pos[u]
                    x+= dx
                    obstruction[y] = max(x+1, obstruction[y])
                    pos[u] = x,y
                    nextlevel += children[u]

                level = nextlevel

        while stack:
            C = stack[-1]
            if len(C) == 0:
                p = stick.pop()
                stack.pop()
                cp = children[p]
                y = o*len(stack)
                if len(cp) == 0:
                    x = obstruction[y]
                    pos[p] = x,y
                else:
                    x = sum([pos[c][0] for c in cp])/(float(len(cp)))
                    pos[p] = x,y
                    ox = obstruction[y]
                    if x < ox:
                        slide(p,ox-x)
                        x = ox
                obstruction[y] = x+1
                continue

            t = C.pop()
            pt = parent[t]

            ct = [u for u in T.neighbors(t) if u != pt]
            for c in ct:
                parent[c] = t
            children[t] = ct

            stack.append([c for c in ct])
            stick.append(t)

        return pos
