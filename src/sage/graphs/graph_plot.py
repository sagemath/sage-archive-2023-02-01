
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
from sage.misc.randstate import current_randstate
from sage.structure.sage_object import SageObject
import sage.graphs.graph_fast as graph_fast
from sage.plot.all import Graphics, scatter_plot, bezier_path, line, arrow, text, circle
from sage.plot.misc import options
from math import sqrt, cos, sin, atan, pi

graphplot_options = {'pos': 'The position dictionary of vertices',
                    'layout': 'A specified layout style-- one of "spring", "circular", "tree".',
                    'vertex_labels': 'Whether or not to draw vertex labels.',
                    'vertex_colors': 'Dictionary of vertex coloring.',
                    'vertex_size': 'The size to draw the vertices.',
                    'vertex_shape': 'The shape to draw the vertices, Currently unavailable for Multi-edged DiGraphs.',
                    'edge_labels': 'Whether or not to draw edge labels.',
                    'edge_style': 'The linestyle of the edges-- one of "solid", "dashed", "dotted", dashdot".',
                    'edge_colors': 'Dictionary of edge coloring.',
                    'color_by_label': 'Whether or not to color the edges by their label values.',
                    'partition': 'A partition of the vertex set.  (Draws each cell of vertices in a different color).',
                    'iterations': 'The number of times to execute the spring layout algorithm.',
                    'loop_size': 'The radius of the smallest loop.',
                    'dist': 'The distance between multiedges.',
                    'max_dist': 'The max distance range to allow multiedges.',
                    'heights': 'Dictionary specifying height (y positions) for vertices.',
                    'tree_root': 'A vertex designation for drawing trees.',
                    'tree_orientation': 'The direction of tree branches-- "up" or "down".',
                    'save_pos': 'Whether or not to save the computed position for the graph.',
                    'graph_border': 'Whether or not to draw a frame around the graph.'}

class GraphPlot(SageObject):
    def __init__(self, graph, options):
        """
        Returns a GraphPlot object, which stores all the parameters needed for
        plotting (Di)Graphs.  A GraphPlot has a plot and show function, as well
        as some functions to set parameters for vertices and edges.  This constructor
        assumes default options are set.  Defaults are shown in the example below.

        EXAMPLE:
            sage: from sage.graphs.graph_plot import GraphPlot
            sage: options = {
            ...     'vertex_size':200,
            ...     'vertex_labels':True,
            ...     'layout':None,
            ...     'edge_style':'solid',
            ...     'edge_colors':'black',
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
        """
        self._plot_components = {}
        self._nodelist = graph.vertices()
        self._graph = graph
        self._options = options
        self.set_pos()
        if self._graph.is_directed() and self._graph.multiple_edges():
            self._multidigraph = True
        else:
            self._multidigraph = False
        self.set_vertices()
        self.set_edges()

    def _repr_(self):
        """
        Returns a string representation of a GraphPlot object.

        EXAMPLE:
            This function is called implicitly by the code below:
            sage: g = Graph({0:[1,2], 2:[3], 4:[0,1]})
            sage: g.graphplot()
            GraphPlot object for Graph on 5 vertices
        """
        return "GraphPlot object for %s"%self._graph

    def set_pos(self):
        """
        Sets the position plotting parameters for this GraphPlot.

        EXAMPLES:
            This function is called implicitly by the code below:
            sage: g = Graph({0:[1,2], 2:[3], 4:[0,1]})
            sage: g.graphplot(save_pos=True, layout='circular')
            GraphPlot object for Graph on 5 vertices

            The following illustrates the format of a position dictionary,
            but due to numerical noise we do not check the values themselves.

            sage: g.get_pos()
            {0: [..., ...],
             1: [..., ...],
             2: [..., ...],
             3: [..., ...],
             4: [..., ...]}

            sage: T = list(graphs.trees(7))
            sage: t = T[3]
            sage: t.plot(heights={0:[0], 1:[4,5,1], 2:[2], 3:[3,6]})
        """
        random = current_randstate().python_random().random

        # Check layouts first to override pos
        if 'layout' in self._options and self._options['layout'] is not None:
            if self._options['layout'] == 'spring':
                self._pos = graph_fast.spring_layout_fast(self._graph, iterations=self._options['iterations'], height=(self._options['heights'] is not None))
            elif self._options['layout'] == 'circular':
                #from math import sin, cos, pi
                n = len(self._nodelist)
                pos = {}
                verts = self._graph.vertices()
                for i in range(n):
                    x = float(cos((pi/2) + ((2*pi)/n)*i))
                    y = float(sin((pi/2) + ((2*pi)/n)*i))
                    pos[verts[i]] = [x,y]
                self._pos = pos
            elif self._options['layout'] == 'tree':
                if not self._graph.is_tree():
                    raise RuntimeError("Cannot use tree layout on this graph: self.is_tree() returns False.")
                n = self._graph.order()
                if not 'tree_root' in self._options:
                    from sage.misc.prandom import randrange
                    root = self._nodelist[randrange(n)]
                else:
                    root = self._options['tree_root']
                # BFS search for heights
                seen = [root]
                queue = [root]
                heights = [-1]*n
                heights[self._nodelist.index(root)] = 0
                while queue:
                    u = queue.pop(0)
                    for v in self._graph.neighbors(u):
                        if v not in seen:
                            seen.append(v)
                            queue.append(v)
                            heights[self._nodelist.index(v)] = heights[self._nodelist.index(u)] + 1
                if self._options['tree_orientation'] == 'down':
                    maxx = max(heights)
                    heights = [maxx-heights[i] for i in xrange(n)]
                heights_dict = {}
                for v in self._nodelist:
                    if not heights_dict.has_key(heights[self._nodelist.index(v)]):
                        heights_dict[heights[self._nodelist.index(v)]] = [v]
                    else:
                        heights_dict[heights[self._nodelist.index(v)]].append(v)
                # Add/overwrite self._options['heights']
                self._options['heights'] = heights_dict

        # No layout-- obtain position dict for all nodes
        elif 'pos' in self._options:
            self._pos = self._options['pos']
        else:
            self._pos = self._graph.get_pos()

        # Set/adjust pos dict for heights
        if self._options['heights'] is not None and self._graph.order() > 0:
            pos = {}
            mmax = max([len(ccc) for ccc in self._options['heights'].values()])
            ymin = min(self._options['heights'].keys())
            ymax = max(self._options['heights'].keys())
            dist = ((ymax-ymin)/(mmax+1.0))
            for height in self._options['heights']:
                num_xs = len(self._options['heights'][height])
                if num_xs == 0:
                    continue
                j = (mmax - num_xs)/2.0
                for k in range(num_xs):
                    pos[self._options['heights'][height][k]] = [dist*(j+k+1) + random()*(dist*0.03), height]
                self._pos = pos

        if not self._pos:
            self._pos = graph_fast.spring_layout_fast(self._graph, iterations=self._options['iterations'], height=(self._options['heights'] is not None))

        # Collect max/min values to add positions if necessary
        if not self._pos:
            self._pos = {}
            xmin = -1
            xmax = 1
            ymin = -1
            ymax = 1
        else:
            xs = [self._pos[v][0] for v in self._pos]
            ys = [self._pos[v][1] for v in self._pos]
            xmin = min(xs)
            xmax = max(xs)
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

        # Check each vertex position is in pos, add position
        # randomly within the plot range if none is defined
        for v in self._nodelist:
            if not v in self._pos:
                self._pos[v] = [xmin + dx*random(), ymin + dy*random()]

        # Save positions to graph if requested in options
        if 'save_pos' in self._options and self._options['save_pos']:
            self._graph.set_pos(self._pos)

    def set_vertices(self, **vertex_options):
        """
        Sets the vertex plotting parameters for this GraphPlot.  This function
        is called by the constructor but can also be called to make updates to
        the vertex options of an existing GraphPlot object.  Note that the
        changes are cumulative.

        EXAMPLES:
            sage: g = Graph({}, loops=True, multiedges=True)
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

        if 'vertex_colors' not in self._options:
            if self._options['partition'] is not None:
                from sage.plot.all import rainbow
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
            if self._multidigraph:
                self._vertex_radius = sqrt(voptions['markersize']/(20500*pi))

        voptions['zorder'] = 7

        if not isinstance(vertex_colors, dict):
            voptions['facecolor'] = vertex_colors
            if self._multidigraph:
                self._plot_components['vertices'] = [circle(center, self._vertex_radius, fill=True, facecolor=vertex_colors) for center in self._pos.values()]
            else:
                self._plot_components['vertices'] = scatter_plot(self._pos.values(), **voptions)
        else:
            # Color list must be ordered:
            pos = []
            colors = []
            for i in vertex_colors:
                pos += [self._pos[j] for j in vertex_colors[i]]
                colors += [i]*len(vertex_colors[i])
            if self._multidigraph:
                self._plot_components['vertices'] = [circle(pos[i], self._vertex_radius, fill=True, facecolor=colors[i]) for i in len(pos)]
            else:
                self._plot_components['vertices'] = scatter_plot(pos, facecolor=colors, **voptions)

        if self._options['vertex_labels']:
            self._plot_components['vertex_labels'] = []
            # TODO: allow text options
            for v in self._nodelist:
                self._plot_components['vertex_labels'].append(text(str(v), self._pos[v], rgbcolor=(0,0,0), zorder=8))

    def set_edges(self, **edge_options):
        """
        Sets the edge (or arrow) plotting parameters for the GraphPlot object.  This
        function is called by the constructor but can also be called to make updates to
        the vertex options of an existing GraphPlot object.  Note that the changes are
        cumulative.

        EXAMPLES:
            sage: g = Graph({}, loops=True, multiedges=True)
            sage: g.add_edges([(0,0,'a'),(0,0,'b'),(0,1,'c'),(0,1,'d'),
            ...     (0,1,'e'),(0,1,'f'),(0,1,'f'),(2,1,'g'),(2,2,'h')])
            sage: GP = g.graphplot(vertex_size=100, edge_labels=True, color_by_label=True, edge_style='dashed')
            sage: GP.set_edges(edge_style='solid')
            sage: GP.plot()
            sage: GP.set_edges(edge_colors='black')
            sage: GP.plot()

            sage: d = DiGraph({}, loops=True, multiedges=True)
            sage: d.add_edges([(0,0,'a'),(0,0,'b'),(0,1,'c'),(0,1,'d'),
            ...     (0,1,'e'),(0,1,'f'),(0,1,'f'),(2,1,'g'),(2,2,'h')])
            sage: GP = d.graphplot(vertex_size=100, edge_labels=True, color_by_label=True, edge_style='dashed')
            sage: GP.set_edges(edge_style='solid')
            sage: GP.plot()
            sage: GP.set_edges(edge_colors='black')
            sage: GP.plot()
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
        else:
            for edge in self._graph.edges(sort=True):
                key = tuple(sorted([edge[0],edge[1]]))
                if key == (edge[0],edge[1]): head = 1
                else: head = 0
                if key in edges_to_draw:
                    edges_to_draw[key].append((edge[2], self._options['edge_colors'], head))
                else:
                    edges_to_draw[key] = [(edge[2], self._options['edge_colors'], head)]

        if edges_to_draw:
            self._plot_components['edges'] = []
        else:
            return

        # Check for multi-edges or loops
        if self._graph.multiple_edges() or self._graph.loops():
            tmp = edges_to_draw.copy()
            dist = self._options['dist']*2
            loop_size = self._options['loop_size']
            max_dist = self._options['max_dist']
            from sage.calculus.calculus import SymbolicVariable
            from sage.calculus.calculus import sqrt as Sqrt
            for (a,b) in tmp:
                if a == b:
                    # Loops
                    distance = dist
                    local_labels = edges_to_draw.pop((a,b))
                    if len(local_labels)*dist > max_dist:
                        distance = max_dist/len(local_labels)
                    curr_loop_size = loop_size
                    for i in range(len(local_labels)):
                        self._plot_components['edges'].append(circle((self._pos[a][0],self._pos[a][1]-curr_loop_size), curr_loop_size, rgbcolor=local_labels[i][1], **eoptions))
                        if labels:
                            self._plot_components['edge_labels'].append(text(local_labels[i][0], (self._pos[a][0], self._pos[a][1]-2*curr_loop_size)))
                        curr_loop_size += distance/4
                elif len(edges_to_draw[(a,b)]) > 1:
                    # Multi-edge
                    local_labels = edges_to_draw.pop((a,b))
                    x = SymbolicVariable('x')
                    d = SymbolicVariable('d')

                    # Compute perpendicular bisector
                    p1 = self._pos[a]
                    p2 = self._pos[b]
                    M = ((p1[0]+p2[0])/2, (p1[1]+p2[1])/2) # midpoint
                    if not p1[1] == p2[1]:
                        S = (p1[0]-p2[0])/(p2[1]-p1[1]) # perp slope
                        y = S*x-S*M[0]+M[1] # perp bisector line

                        # f,g are functions of distance d to determine x values
                        # on line y at d from point M
                        f = Sqrt(d**2/(1+S**2)) + M[0]
                        g = -Sqrt(d**2/(1+S**2)) + M[0]

                        odd_x = f
                        even_x = g
                        if p1[0] == p2[0]:
                            odd_y = lambda d : M[1]
                            even_y = odd_y
                        else:
                            odd_y = y(f)
                            even_y = y(g)
                    else:
                        odd_x = lambda d : M[0]
                        even_x = odd_x
                        odd_y = M[1] + d
                        even_y = M[1] - d

                    # We now have the control points for each bezier curve
                    # in terms of distance parameter d.
                    # Also note that the label for each edge should be drawn at d/2.
                    # (This is because we're using the perp bisectors).
                    distance = dist
                    if len(local_labels)*dist > max_dist:
                        distance = max_dist/len(local_labels)
                    for i in range(len(local_labels)/2):
                        k = (i+1.0)*distance
                        if self._multidigraph:
                            C,D = self._polar_hack_for_multidigraph(p1, [odd_x(k),odd_y(k)], self._vertex_radius)
                            E,F = self._polar_hack_for_multidigraph([odd_x(k),odd_y(k)], p2, self._vertex_radius)
                            self._plot_components['edges'].append(arrow(path=[[C,[odd_x(k),odd_y(k)],F]], head=local_labels[2*i][2], zorder=1, rgbcolor=local_labels[2*i][1], **eoptions))
                            self._plot_components['edges'].append(arrow(path=[[C,[even_x(k),even_y(k)],F]], head=local_labels[2*i+1][2], zorder=1, rgbcolor=local_labels[2*i+1][1], **eoptions))
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
            if self._multidigraph:
                C,D = self._polar_hack_for_multidigraph(self._pos[a], self._pos[b], self._vertex_radius)
                self._plot_components['edges'].append(arrow(C,D, rgbcolor=edges_to_draw[(a,b)][0][1], head=edges_to_draw[(a,b)][0][2], **eoptions))
                if labels:
                    self._plot_components['edge_labels'].append(text(str(edges_to_draw[(a,b)][0][0]),[(C[0]+D[0])/2, (C[1]+D[1])/2]))
            elif dir:
                self._plot_components['edges'].append(arrow(self._pos[a],self._pos[b], rgbcolor=edges_to_draw[(a,b)][0][1], arrowshorten=self._arrowshorten, head=edges_to_draw[(a,b)][0][2], **eoptions))
            else:
                self._plot_components['edges'].append(line([self._pos[a],self._pos[b]], rgbcolor=edges_to_draw[(a,b)][0][1], **eoptions))
            if labels and not self._multidigraph:
                self._plot_components['edge_labels'].append(text(str(edges_to_draw[(a,b)][0][0]),[(self._pos[a][0]+self._pos[b][0])/2, (self._pos[a][1]+self._pos[b][1])/2]))

    def _polar_hack_for_multidigraph(self, A, B, VR):
        """
        Helper function to quickly compute the two points of intersection of a line
        segment from A to B (xy tuples) and circles centered at A anb B, both with
        radius VR.  Returns a tuple of xy tuples representing the two points.

        EXAMPLE:
            sage: d = DiGraph({}, loops=True, multiedges=True)
            sage: d.add_edges([(0,0,'a'),(0,0,'b'),(0,1,'c'),(0,1,'d'),
            ...     (0,1,'e'),(0,1,'f'),(0,1,'f'),(2,1,'g'),(2,2,'h')])
            sage: GP = d.graphplot(vertex_size=100, edge_labels=True, color_by_label=True, edge_style='dashed')
            sage: GP._polar_hack_for_multidigraph((0,1),(1,1),.1)
            ([0.10..., 1.00...], [0.90..., 1.00...])
        """
        D = [B[i]-A[i] for i in range(2)]
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

        EXAMPLE:
            sage: C = graphs.CubeGraph(8)
            sage: P = C.graphplot(vertex_labels=False, vertex_size=0, graph_border=True)
            sage: P.show()
        """
        self.plot().show(**kwds)

    def plot(self, **kwds):
        """
        Returns a graphics object representing the (di)graph.

        INPUT:
            pos -- an optional positioning dictionary
            layout -- what kind of layout to use, takes precedence over pos
                'circular' -- plots the graph with vertices evenly distributed
                    on a circle
                'spring' -- uses the traditional spring layout, using the
                    graph's current positions as initial positions
                'tree' -- the (di)graph must be a tree. One can specify the root
                    of the tree using the keyword tree_root, otherwise a root
                    will be selected at random. Then the tree will be plotted in
                    levels, depending on minimum distance for the root.
            vertex_labels -- whether to print vertex labels
            edge_labels -- whether to print edge labels. By default, False,
                but if True, the result of str(l) is printed on the edge for
                each label l. Labels equal to None are not printed (to set edge
                labels, see set_edge_label).
            vertex_size -- size of vertices displayed
            vertex_shape -- the shape to draw the vertices (Not available for
                multiedge digraphs.
            graph_border -- whether to include a box around the graph
            vertex_colors -- optional dictionary to specify vertex colors: each
                key is a color recognizable by matplotlib, and each corresponding
                entry is a list of vertices. If a vertex is not listed, it looks
                invisible on the resulting plot (it doesn't get drawn).
            edge_colors -- a dictionary specifying edge colors: each key is a
                color recognized by matplotlib, and each entry is a list of edges.
            partition -- a partition of the vertex set. if specified, plot will
                show each cell in a different color. vertex_colors takes precedence.
            scaling_term -- default is 0.05. if vertices are getting chopped off,
                increase; if graph is too small, decrease. should be positive, but
                values much bigger than 1/8 won't be useful unless the vertices
                are huge
            talk -- if true, prints large vertices with white backgrounds so that
                labels are legible on slides
            iterations -- how many iterations of the spring layout algorithm to
                go through, if applicable
            color_by_label -- if True, color edges by their labels
            heights -- if specified, this is a dictionary from a set of
                floating point heights to a set of vertices
            edge_style -- keyword arguments passed into the
                edge-drawing routine.  This currently only works for
                directed graphs, since we pass off the undirected graph to
                networkx
            tree_root -- a vertex of the tree to be used as the root for
                the layout="tree" option. If no root is specified, then one
                is chosen at random. Ignored unless layout='tree'.
            tree_orientation -- "up" or "down" (default is "down").
                If "up" (resp., "down"), then the root of the tree will
                appear on the bottom (resp., top) and the tree will grow
                upwards (resp. downwards). Ignored unless layout='tree'.
            save_pos -- save position computed during plotting

        EXAMPLES:
            sage: from sage.graphs.graph_plot import graphplot_options
            sage: list(sorted(graphplot_options.iteritems()))
            [('color_by_label', 'Whether or not to color the edges by their label values.'),
            ('dist', 'The distance between multiedges.'),
            ('edge_colors', 'Dictionary of edge coloring.'),
            ('edge_labels', 'Whether or not to draw edge labels.'),
            ('edge_style', 'The linestyle of the edges-- one of "solid", "dashed", "dotted", dashdot".'),
            ('graph_border', 'Whether or not to draw a frame around the graph.'),
            ('heights', 'Dictionary specifying height (y positions) for vertices.'),
            ('iterations', 'The number of times to execute the spring layout algorithm.'),
            ('layout', 'A specified layout style-- one of "spring", "circular", "tree".'),
            ('loop_size', 'The radius of the smallest loop.'),
            ('max_dist', 'The max distance range to allow multiedges.'),
            ('partition', 'A partition of the vertex set.  (Draws each cell of vertices in a different color).'),
            ('pos', 'The position dictionary of vertices'),
            ('save_pos', 'Whether or not to save the computed position for the graph.'),
            ('tree_orientation', 'The direction of tree branches-- "up" or "down".'),
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

            sage: G = graphs.HeawoodGraph()
            sage: for u,v,l in G.edges():
            ...    G.set_edge_label(u,v,'(' + str(u) + ',' + str(v) + ')')
            sage: G.graphplot(edge_labels=True).show()

            sage: D = DiGraph( { 0: [1, 10, 19], 1: [8, 2], 2: [3, 6], 3: [19, 4], 4: [17, 5], 5: [6, 15], 6: [7], 7: [8, 14], 8: [9], 9: [10, 13], 10: [11], 11: [12, 18], 12: [16, 13], 13: [14], 14: [15], 15: [16], 16: [17], 17: [18], 18: [19], 19: []}, implementation='networkx' )
            sage: for u,v,l in D.edges():
            ...    D.set_edge_label(u,v,'(' + str(u) + ',' + str(v) + ')')
            sage: D.graphplot(edge_labels=True, layout='circular').show()

            sage: from sage.plot.plot import rainbow
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

            sage: D = DiGraph(multiedges=True)
            sage: for i in range(5):
            ...     D.add_edge((i,i+1,'a'))
            ...     D.add_edge((i,i-1,'b'))
            sage: D.graphplot(edge_labels=True,edge_colors=D._color_by_label()).plot()

            sage: g = Graph({}, loops=True, multiedges=True)
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
        return G


