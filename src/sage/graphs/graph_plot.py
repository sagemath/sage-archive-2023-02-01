
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

graphplot_options = {'pos': 'The position dictionary of vertices',
                    'layout': 'A specified layout style-- one of "spring", "circular", "tree".',
                    'vertex_labels': 'Whether or not to draw vertex labels.',
                    'vertex_colors': 'Dictionary of vertex coloring.',
                    'vertex_size': 'The size to draw the vertices.',
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
        TODO: docstring
        """

        self._plot_components = {}
        self._nodelist = graph.vertices()
        self._graph = graph
        self._options = options
        self.set_pos()
        self.set_vertices()

        if not self._graph.is_directed():
            self.set_edges()
        else:
            self.set_arrows()

    def _repr_(self):
        return "GraphPlot object for %s"%self._graph

    def set_pos(self):
        """
        TODO
        """
        random = current_randstate().python_random().random

        # Check layouts first to override pos
        if 'layout' in self._options and self._options['layout'] is not None:
            if self._options['layout'] == 'spring':
                self._pos = graph_fast.spring_layout_fast(self._graph, iterations=self._options['iterations'], height=(self._options['heights'] is not None))
            elif self._options['layout'] == 'circular':
                from math import sin, cos, pi
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
        TODO
        """
        # Handle base vertex options
        voptions = vertex_options

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

        voptions['zorder'] = 7

        # TODO: allow advanced (**vertex_options) options to be added by user

        if not isinstance(vertex_colors, dict):
            voptions['facecolor'] = vertex_colors
            self._plot_components['vertices'] = scatter_plot(self._pos.values(), **voptions)
        else:
            # Color list must be ordered:
            pos = []
            colors = []
            for i in vertex_colors:
                pos += [self._pos[j] for j in vertex_colors[i]]
                colors += [i]*len(vertex_colors[i])
            self._plot_components['vertices'] = scatter_plot(pos, facecolor=colors, **voptions)

        if self._options['vertex_labels']:
            self._plot_components['vertex_labels'] = []
            # TODO: allow text options
            for v in self._nodelist:
                self._plot_components['vertex_labels'].append(text(str(v), self._pos[v], rgbcolor=(0,0,0), zorder=8))

    def set_edges(self, **edge_options):
        """
        TODO: will eventually want to allow user to entirely specify edge components
        ie.: give actual curve parameters for each edge
        """
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
                    if (edge[0],edge[1]) in edges_to_draw:
                        edges_to_draw[(edge[0],edge[1])].append((edge[2], color))
                    else:
                        edges_to_draw[(edge[0],edge[1])] = [(edge[2], color)]
        else:
            for edge in self._graph.edges(sort=True):
                if (edge[0],edge[1]) in edges_to_draw:
                    edges_to_draw[(edge[0],edge[1])].append((edge[2], self._options['edge_colors']))
                else:
                    edges_to_draw[(edge[0],edge[1])] = [(edge[2], self._options['edge_colors'])]

        if edges_to_draw:
            self._plot_components['edges'] = []

        # Check for multi-edges or loops
        if self._graph.multiple_edges() or self._graph.loops():
            tmp = edges_to_draw.copy()
            dist = self._options['dist']*2
            loop_size = self._options['loop_size']
            max_dist = self._options['max_dist']
            from sage.calculus.calculus import SymbolicVariable, sqrt
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
                        # from math import sqrt
                        f = sqrt(d**2/(1+S**2)) + M[0]
                        g = -sqrt(d**2/(1+S**2)) + M[0]

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
                        self._plot_components['edges'].append(bezier_path([[p1,[odd_x(k),odd_y(k)],p2]],zorder=1, rgbcolor=local_labels[2*i][1], **eoptions))
                        self._plot_components['edges'].append(bezier_path([[p1,[even_x(k),even_y(k)],p2]],zorder=1, rgbcolor=local_labels[2*i+1][1], **eoptions))
                        if labels:
                            j = k/2.0
                            self._plot_components['edge_labels'].append(text(local_labels[2*i][0],[odd_x(j),odd_y(j)]))
                            self._plot_components['edge_labels'].append(text(local_labels[2*i+1][0],[even_x(j),even_y(j)]))
                    if len(local_labels)%2 == 1:
                        edges_to_draw[(a,b)] = [local_labels[-1]] # draw line for last odd

        # Add lines for remainder of edges
        for (a,b) in edges_to_draw:
            self._plot_components['edges'].append(line([self._pos[a],self._pos[b]], rgbcolor=edges_to_draw[(a,b)][0][1], **eoptions))
            if labels:
                self._plot_components['edge_labels'].append(text(str(edges_to_draw[(a,b)][0][0]),[(self._pos[a][0]+self._pos[b][0])/2, (self._pos[a][1]+self._pos[b][1])/2]))

    def set_arrows(self, **arrow_options):
        pass

    def show(self, **kwds):
        self.plot().show(**kwds)

    def plot(self, **kwds):
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


