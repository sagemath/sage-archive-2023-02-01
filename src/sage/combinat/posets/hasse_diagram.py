r"""
Hasse diagrams of posets
"""
#*****************************************************************************
#       Copyright (C) 2008 Peter Jipsen <jipsen@chapman.edu>,
#                          Franco Saliola <saliola@gmail.com>
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

import __builtin__
from sage.structure.parent_base import ParentWithBase
from sage.structure.element import Element
from sage.rings.integer import Integer
from sage.graphs.graph import Graph, DiGraph
from sage.misc.sage_eval import sage_eval
from sage.matrix.constructor import matrix
from sage.rings.finite_field import FiniteField
from sage.rings.integer_ring import IntegerRing, ZZ
from sage.matrix.constructor import matrix
from sage.misc.misc import uniq

class HasseDiagram(DiGraph):
    """
    The Hasse diagram of a poset. This is just a transitively-reduced,
    directed, acyclic graph without loops or multiple edges.

    NOTE: We assume that range(n) is a linear extension of the poset.

    That is, range(n) is the vertex set and a topological sort of the
    digraph.

    This should not be called directly, use Poset instead.
    """
    def __repr__(self):
        return "Hasse diagram of a poset containing %s elements"%self.order()

    def linear_extension(self):
        # Recall: we assume range(n) is a linear extension.
        return range(len(self))

    def linear_extensions(self):
        return self.topological_sort_generator()

    def is_linear_extension(self,lin_ext=None):
        if lin_ext is None or lin_ext == range(len(self)):
            for x,y in self.cover_relations_iterator():
                if not x < y:
                    return False
            return True
        else:
            for x,y in self.cover_relations_iterator():
                if not lin_ext.index(x) < lin_ext.index(y):
                    return False
            return True

    def level_sets(self):
        """
        Returns a list l such that l[i+1] is the set of minimal
        elements of the poset obtained by removing the elements
        in l[0], l[1], ..., l[i].

        EXAMPLES:
            sage: P = Poset({0:[1,2],1:[3],2:[3],3:[]})
            sage: [len(x) for x in P.level_sets()]
            [1, 2, 1]

            sage: Q = Poset({0:[1,2], 1:[3], 2:[4], 3:[4]})
            sage: [len(x) for x in Q.level_sets()]
            [1, 2, 1, 1]
        """
        G = self.copy()
        Levels = []
        while G.vertices() != []:
            indegs = G.in_degree(labels=True)
            new_level = [x for x in indegs if indegs[x]==0]
            Levels.append(new_level)
            G.delete_vertices(new_level)
        return Levels

    def plot(self, label_elements=True, element_labels=None,
            label_font_size=12,label_font_color='black',
            vertex_size=300, vertex_colors=None):
        """
	Returns a Graphics object corresponding to the Hasse diagram.

        EXAMPLES:
            sage: uc = [[2,3], [], [1], [1], [1], [3,4]]
            sage: elm_lbls = Permutations(3).list()
            sage: P = Poset(uc,elm_lbls)
            sage: H = P.hasse_diagram()
            sage: levels = H.level_sets()
            sage: heights = dict([[i, levels[i]] for i in range(len(levels))])
            sage: type(H.plot(label_elements=True))
            <class 'sage.plot.plot.Graphics'>
        """
        # If we have no elements: plot an empty poset (blank plot).
        if self.order() == 0:
            return networkx_plot(self.networkx_graph())

        # Set element_labels to default to the vertex set.
        if element_labels is None:
            element_labels = dict(zip(self.vertices(),self.vertices()))

        # Compute dictionary of heights.
        if self._pos is None:
            # Determine heights of vertices based on levels.
            levels = self.level_sets()
            def sort_fcn(x,y):
                if element_labels[x] > element_labels[y]:
                    return 1
                elif element_labels[x] < element_labels[y]:
                    return -1
                else:
                    return 0
            levels = [sorted(z,sort_fcn) for z in levels]
            heights = dict([[i, levels[i]] for i in range(len(levels))])


            # Set positions for vertices based on height.
            # (Code cripped from graph.plot().)
            pos = {}
            mmax = max([len(ccc) for ccc in heights.values()])
            ymin = min(heights.keys())
            ymax = max(heights.keys())
            if ymin==ymax:
                dist = 0.5
            else:
                dist = ((ymax-ymin)/(mmax+1.0))
            for height in heights:
                num_xes = len(heights[height])
                if num_xes == 0: continue
                j = (mmax - num_xes)/2.0
                for k in range(num_xes):
                    pos[heights[height][k]] = [dist*(j+k+1),height]
        else:
            pos = self._pos

        # Create the underlying graph.
        graph = self.to_undirected()

        # Set default vertex colors.
        if vertex_colors is None:
            vertex_colors = {'#FFFFFF':graph.vertices()}

        # Create the networkx plot.
        from sage.plot.plot import networkx_plot
        G = networkx_plot(graph.networkx_graph(), pos=pos, \
                vertex_size=vertex_size, vertex_labels=False, \
                vertex_colors=vertex_colors)

        # Add vertex labels to the graph.
        if label_elements:
            from sage.plot.plot import text, polygon, Graphics
            K = Graphics()
            for v in graph.vertices():
                K += text(str(element_labels[v]), pos[v],
                        rgbcolor=label_font_color, fontsize=label_font_size)
            K.axes_range(xmin=G.xmin(), xmax=G.xmax(), ymin=G.ymin(), ymax=G.ymax())
            G += K
        G.axes(False)

        # Return the graph
        return G

    def show(self, label_elements=True, element_labels=None,
            label_font_size=12,label_font_color='black',
            vertex_size=300, vertex_colors=None,**kwds):
        """
	Shows the Graphics object corresponding to the Hasse diagram.
	Optionally, it is labelled.

        INPUT:
            label_elements -- whether to display element labels
            element_labels -- a dictionary of element labels

        EXAMPLES:
            sage: uc = [[2,3], [], [1], [1], [1], [3,4]]
            sage: elm_lbls = Permutations(3).list()
            sage: P = Poset(uc,elm_lbls)
            sage: H = P.hasse_diagram()
            sage: levels = H.level_sets()
            sage: heights = dict([[i, levels[i]] for i in range(len(levels))])
            sage: H.show(label_elements=True)
        """
        self.plot(label_elements=label_elements, element_labels=element_labels,
            label_font_size=label_font_size,label_font_color=label_font_color,
            vertex_size=vertex_size, vertex_colors=vertex_colors).show(**kwds)

    def cover_relations_iterator(self):
        for u,v,l in self.edge_iterator():
            yield (u,v)

    def cover_relations(self,element=None):
        return [c for c in self.cover_relations_iterator()]

    def is_lequal(self, i, j):
        """
        Returns True if i is less than or equal to j in the poset, and False
        otherwise.

        EXAMPLES:
        """
        return i == j or \
                (i < j and j in self.breadth_first_search(i))

    def is_less_than(self, x, y):
        if x == y:
            return False
        else:
            return self.is_lequal(x,y)

    def is_gequal(self, x, y):
        return self.is_lequal(y,x)

    def is_greater_than(self, x, y):
        return self.is_lequal(y,x)

    def minimal_elements(self):
        """
        Returns a list of the minimal elements of the poset.

        EXAMPLES:
            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P(0) in P.minimal_elements()
            True
            sage: P(1) in P.minimal_elements()
            True
            sage: P(2) in P.minimal_elements()
            True
        """
        indegs = self.in_degree(labels=True)
        return [x for x in indegs if indegs[x]==0]

    def maximal_elements(self):
        """
        Returns a list of the maximal elements of the poset.

        EXAMPLES:
            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P.maximal_elements()
            [4]
        """
        outdegs = self.out_degree(labels=True)
        return [x for x in outdegs if outdegs[x]==0]

    def bottom(self):
        """
        Returns the bottom element of the poset, if it exists.

        EXAMPLES:
            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P.bottom() is None
            True
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.bottom()
            0
        """
        min_elms = self.minimal_elements()
        if len(min_elms) == 1: return min_elms[0]
        return None

    def has_bottom(self):
        """
        Returns True if the poset has a unique minimal element.

        EXAMPLES:
            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P.has_bottom()
            False
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.has_bottom()
            True
        """
        if self.bottom() is not None: return True
        return False

    def top(self):
        """
        Returns the top element of the poset, if it exists.

        EXAMPLES:
            sage: P = Poset({0:[3],1:[3],2:[3],3:[4,5],4:[],5:[]})
            sage: P.top() is None
            True
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.top()
            1
        """
        max_elms = self.maximal_elements()
        if len(max_elms) == 1: return max_elms[0]
        return None

    def has_top(self):
        """
        Returns True if the poset contains a unique maximal element,
        and False otherwise.

        EXAMPLES:
            sage: P = Poset({0:[3],1:[3],2:[3],3:[4,5],4:[],5:[]})
            sage: P.has_top()
            False
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.has_top()
            True
        """
        if not self.top() is None: return True
        return False

    def is_bounded(self):
        """
        Returns True if the poset contains a unique maximal element
        and a unique minimal element, and False otherwise.

        EXAMPLES:
            sage: P = Poset({0:[3],1:[3],2:[3],3:[4,5],4:[],5:[]})
            sage: P.is_bounded()
            False
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.is_bounded()
            True
        """
        return self.has_top() and self.has_bottom()

    def is_chain(self):
        """
        Returns True if the poset is totally ordered, and False
        otherwise.

        EXAMPLES:
            sage: L = Poset({0:[1],1:[2],2:[3],3:[4]})
            sage: L.is_chain()
            True
            sage: V = Poset({0:[1,2]})
            sage: V.is_chain()
            False
        """
        outdegs = self.out_degree()
        outdegs.remove(0)
        if len(set(outdegs))==1: return True
        return False

    def dual(self):
        """
        Returns a poset that is dual to the given poset.

        EXAMPLE:
            sage: P = Poset([[1,2],[4],[3],[4],[]])
            sage: P.dual()
            Finite poset containing 5 elements
        """
        dual_graph = self.reverse()
        dual_lin_ext = copy(self.linear_extension())
        dual_lin_ext = dual_lin_ext.reverse()
        return HasseDiagram(dual_graph,dual_lin_ext)

    def interval(self, x, y):
        """
        Returns a list of the elements z such that x <= z <= y. The
        order is that induced by the ordering in self.linear_extension.

        INPUT:
            x -- any element of the poset
            y -- any element of the poset

        EXAMPLES:
            sage: uc = [[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]]
            sage: dag = DiGraph(dict(zip(range(len(uc)),uc)))
            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram(dag)
            sage: I = set([2,5,6,4,7])
            sage: I == set(H.interval(2,7))
            True
        """
        return [z for z in range(self.order())[x:y+1] if
                self.is_lequal(x,z) and self.is_lequal(z,y)]

    def closed_interval(self, x, y):
        """
        Returns a list of the elements z such that x <= z <= y. The
        order is that induced by the ordering in self.linear_extension.

        EXAMPLES:
            sage: uc = [[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]]
            sage: dag = DiGraph(dict(zip(range(len(uc)),uc)))
            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram(dag)
            sage: set([2,5,6,4,7]) == set(H.closed_interval(2,7))
            True
        """
        return self.interval(x,y)

    def open_interval(self, x, y):
        """
        Returns a list of the elements z such that x < z < y. The
        order is that induced by the ordering in self.linear_extension.

        EXAMPLES:
            sage: uc = [[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]]
            sage: dag = DiGraph(dict(zip(range(len(uc)),uc)))
            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram(dag)
            sage: set([5,6,4]) == set(H.open_interval(2,7))
            True
            sage: H.open_interval(7,2)
            []
        """
        ci = self.interval(x,y)
        if len(ci) == 0:
            return []
        else:
            return ci[1:-1]

    def rank_function(self):
        r"""
        Returns a rank function of the poset, if it exists.

        A \textit{rank function} of a poset $P$ is a function $r$ from
        that maps elements of $P$ to integers and satisfies: $r(x) =
        r(y) + 1$ if $x$ covers $y$.

        EXAMPLES:
            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: P.rank_function() is not None
            True
            sage: r = P.rank_function()
            sage: for u,v in P.cover_relations_iterator():
            ...    if r(v) != r(u) + 1:
            ...        print "Bug in rank_function!"

            sage: Q = Poset([[1,2],[4],[3],[4],[]])
            sage: Q.rank_function() is None
            True
        """
        if hasattr(self,"_rank_function"):
            return self._rank_function
        levels = self.level_sets()
        rank_fcn = {}
        for i in range(len(levels)):
            for x in levels[i]:
                rank_fcn[x]=i
        for e in self.cover_relations_iterator():
            if rank_fcn[e[1]]-rank_fcn[e[0]] != 1:
                return None
        self._rank_function = lambda u: rank_fcn[u]
        return self._rank_function

    def rank(self,element=None):
        if element is None:
            return len(self.level_sets())-1
        else:
            return self.rank_function()(element)

    def is_ranked(self):
        r"""
        Returns True if the poset is ranked, and False otherwise.

        A poset is \texit{ranked} if it admits a rank function.

        EXAMPLES:
            sage: P = Poset([[1],[2],[3],[4],[]])
            sage: P.is_ranked()
            True
            sage: Q = Poset([[1,5],[2,6],[3],[4],[],[6,3],[4]])
            sage: Q.is_ranked()
            False
        """
        return bool(self.rank_function())

    def is_graded(self):
        r"""
        Returns True if the poset is graded, and False otherwise.

        A poset is \texit{graded} if it admits a rank function.

        EXAMPLES:
            sage: P = Poset([[1],[2],[3],[4],[]])
            sage: P.is_graded()
            True
            sage: Q = Poset([[1,5],[2,6],[3],[4],[],[6,3],[4]])
            sage: Q.is_graded()
            False
        """
        return self.is_ranked()

    def covers(self,x,y):
        """
        Returns True if y covers x and False otherwise.

        EXAMPLES:
            sage: Q = Poset([[1,5],[2,6],[3],[4],[],[6,3],[4]])
            sage: Q.covers(Q(1),Q(6))
            True
            sage: Q.covers(Q(1),Q(4))
            False
        """
        return self.has_edge(x,y)

    def upper_covers_iterator(self,element):
        for x in self.successor_iterator(element):
            yield x

    def lower_covers_iterator(self,element):
        for x in self.predecessor_iterator(element):
            yield x

    def size(self):
        """
        Returns the number of elements in the poset.

        EXAMPLES:
            sage: Poset([[1,2,3],[4],[4],[4],[]]).size()
            5
        """
        return self.order()

    def mobius_function(self,i,j): # dumb algorithm
        r"""
        Returns the value of the M\"obius function of the poset
        on the elements x and y.

        EXAMPLES:
            sage: P = Poset([[1,2,3],[4],[4],[4],[]])
            sage: H = P.hasse_diagram()
            sage: H.mobius_function(0,4)
            2
            sage: for u,v in P.cover_relations_iterator():
            ...    if P.mobius_function(u,v) != -1:
            ...        print "Bug in mobius_function!"
        """
        try:
            return self._mobius_function_values[(i,j)]
        except AttributeError:
            self._mobius_function_values = {}
            return self.mobius_function(i,j)
        except KeyError:
            if i == j:
                self._mobius_function_values[(i,j)] = 1
            elif i > j:
                self._mobius_function_values[(i,j)] = 0
            else:
                ci = self.closed_interval(i,j)
                if len(ci) == 0:
                    self._mobius_function_values[(i,j)] = 0
                else:
                    self._mobius_function_values[(i,j)] = \
                     -sum([self.mobius_function(i,k) for k in ci[:-1]])
        return self._mobius_function_values[(i,j)]

    def mobius_function_matrix(self):
        if not hasattr(self,'_mobius_function_matrix'):
            self._mobius_function_matrix = self.lequal_matrix(ring=ZZ).inverse()

        # Redefine self.mobius_function
        def mobius_function(i,j):
            return self._mobius_function_matrix[i,j]
        self.mobius_function = mobius_function

        return self._mobius_function_matrix

    def order_filter(self,elements):
        """
        Returns the order filter generated by a list of elements.

        I is an order filter if it satisfies:
          x in I  and y >= x implies y in I

        EXAMPLES:
            sage: H = BooleanLattice(4).hasse_diagram()
            sage: H.order_filter([3,8])
            [3, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        """
        of = []
        for i in elements:
            for j in self.breadth_first_search(i):
                of.append(j)
        return uniq(of)

    def principal_order_filter(self, i):
        """
        Returns the order filter generated by i.

        EXAMPLES:
            sage: H = BooleanLattice(4).hasse_diagram()
            sage: H.principal_order_filter(2)
            [2, 3, 6, 7, 10, 11, 14, 15]
        """
        return self.order_filter([i])

    def order_ideal(self,elements):
        """
        Returns the order ideal generated by a list of elements.

        I is an order ideal if it satisfies:
          x in I  and y <= x implies y in I

        EXAMPLES:
            sage: H = BooleanLattice(4).hasse_diagram()
            sage: H.order_ideal([7,10])
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 10]
        """
        H = self.copy().reverse()
        oi = []
        for i in elements:
            for j in H.breadth_first_search(i):
                oi.append(j)
        return uniq(oi)

    def principal_order_ideal(self, i):
        """
        Returns the order ideal generated by i.

        EXAMPLES:
            sage: H = BooleanLattice(4).hasse_diagram()
            sage: H.principal_order_ideal(6)
            [0, 2, 4, 6]
        """
        return self.order_ideal([i])

    def lequal_matrix(self,ring=ZZ,sparse=True):
        """
        Computes a matrix whose [i,j] entry is 1 if
        self.linear_extension[i] < self.linear_extension[j] and 0
        otherwise. The matrix is stored in the attribute _leq_matrix
        and the method __lt__ is redefined to use this matrix.

        EXAMPLES:
            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: H = P.hasse_diagram()
            sage: hasattr(H,'_leq_matrix')
            False
            sage: H.lequal_matrix()
            [1 1 1 1 1 1 1 1]
            [0 1 0 1 0 0 0 1]
            [0 0 1 1 1 0 1 1]
            [0 0 0 1 0 0 0 1]
            [0 0 0 0 1 0 0 1]
            [0 0 0 0 0 1 1 1]
            [0 0 0 0 0 0 1 1]
            [0 0 0 0 0 0 0 1]
            sage: hasattr(H,'_leq_matrix')
            True
        """
        # If we've already computed the matrix, then return it.
        if hasattr(self,"_leq_matrix") and self._leq_matrix.base_ring() == ring:
            return self._leq_matrix

        # Create the matrix and store it as self._leq_matrix.
        n = self.order()
        D = {}
        for i in range(n):
            for v in self.breadth_first_search(i):
                D[(i,v)] = 1
        self._leq_matrix = matrix(ring, n, n, D, sparse=sparse)

        # Redefine self.is_lequal
        def is_lequal(i, j):
            """
            Returns True if i is less than or equal to j in the poset, and False
            otherwise.

            EXAMPLES:
            """
            return bool(self._leq_matrix[i,j])
        self.is_lequal=is_lequal

        # Return the matrix.
        return self._leq_matrix

    def meet_matrix(self): # modelled after Freese-Jezek-Nation p217
        if hasattr(self,'_meet'):
            return self._meet
        n = self.size()
        meet = [[0 for x in range(n)] for x in range(n)]
        le = self.lequal_matrix()
        for i in range(n): le[i,i] = 1
        if not all([le[0,x]==1 for x in range(n)]):
            raise ValueError, "Not a meet-semilattice: no bottom element."
        lc = [[y[0] for y in self.incoming_edges([x])] for x in range(n)]
        S = []
        for x in range(n): # x=x_k
            meet[x][x] = x
            for y in S:
                T = []
                for z in lc[x]:
                    T.append(meet[y][z]) # T = {x_i \wedge z : z>-x_k}
                q = T[0]
                for z in T[1:]:
                    if z>q: q = z
                for z in T:
                    if not le[z,q]:
                        raise ValueError, "No meet for x=%s y=%s"%(x,y)
                meet[x][y] = q
                meet[y][x] = q
            S.append(x)
        self._meet = matrix(ZZ,meet)
        return self._meet

    def is_meet_semilattice(self):
        try:
            self.meet_matrix()
        except ValueError:
            return False
        else:
            return True

    def join_matrix(self): # Freese-Jezek-Nation p217
        if hasattr(self,'_join'): return self._join
        n = self.size()
        join = [[0 for x in range(n)] for x in range(n)]
        le = self.lequal_matrix()
        for i in range(n): le[i,i] = 1
        if not all([le[x][n-1]==1 for x in range(n)]):
            raise ValueError, "Not a join-semilattice: no top element."
        uc = [sorted([n-1-y[1] for y in self.outgoing_edges([x])]) for
                x in reversed(range(n))]
        S = []
        for x in range(n): # x=x_k
            join[x][x] = x
            for y in S:
                T = []
                for z in uc[x]:
                    T.append(join[y][z]) # T = {x_i \vee z : z>-x_k}
                q = T[0]
                for z in T[1:]:
                    if z>q: q = z
                for z in T:
                    if not le[n-1-q][n-1-z]:
                        raise ValueError, "No join for x=%s y=%s"%(x,y)
                join[x][y] = q
                join[y][x] = q
            S.append(x)
        self._join = matrix(ZZ,[[n-1-join[n-1-x][n-1-y] for y in range(n)] for x in range(n)])
        return self._join

    def is_join_semilattice(self):
        try:
            self.join_matrix()
        except ValueError:
            return False
        else:
            return True

    def is_distributive_lattice(self): # a dumb algorithm...
        try:
            join = self.join_matrix()
            meet = self.meet_matrix()
        except ValueError:
            return False
        n = self.size()
        for x in range(n):
            for y in range(n):
                for z in range(n):
                    if meet[x,join[y,z]] != join[meet[x,y],meet[x,z]]:
                        return False
        return True

    def is_distributive_lattice_fastest(self): # still a dumb algorithm...
        try:
            jn = self.join_matrix()
            mt = self.meet_matrix()
        except ValueError:
            return False
        n = len(jn)
        for x in range(n):
            for y in range(n):
                for z in range(n):
                    if mt[x][jn[y][z]]!=jn[mt[x][y]][mt[x][z]]: return False
        return True

    def is_complemented_lattice(self):
        try:
            jn = self.join_matrix()
            mt = self.meet_matrix()
        except ValueError:
            return False
        n = self.size()
        c = [-1 for x in range(n)]
        for x in range(n):
            for y in range(x,n):
                if jn[x][y]==n-1 and mt[x][y]==0:
                    c[x]=y
                    c[y]=x
        return all([c[x]!=-1 for x in range(n)])

    def complements(self):
        jn = self.join_matrix()
        mt = self.meet_matrix()
        n = self.size()
        c = [-1 for x in range(n)]
        for x in range(n):
            for y in range(x,n):
                if jn[x][y]==n-1 and mt[x][y]==0:
                    c[x]=y
                    c[y]=x
        return c

    def antichains(self): # Freese-Jezek-Nation p226
	"""
        Returns a list of all antichains of the poset.

        An antichain of a poset is a collection of elements of the poset
        that are pairwise incomparable.

        EXAMPLES:
            sage: PentagonPoset().hasse_diagram().antichains()
            [[], [0], [1], [2], [3], [4], [1, 2], [1, 3]]
            sage: AntichainPoset(3).hasse_diagram().antichains()
	    [[], [0], [1], [2], [1, 2], [0, 1], [0, 2], [0, 1, 2]]
	    sage: ChainPoset(3).hasse_diagram().antichains()
	    [[], [0], [1], [2]]
	"""
        S = [[]]
        leq = self.lequal_matrix()
        def antichains_rec(A,x,T,S):
            # A is an antichain. Each element of A is incomparable with
            # x and with each element of T. Also x is not in T.
            Ap = A+[x]
            S.append(Ap)
            if T!=[]:
                antichains_rec(A,T[0],T[1:],S)
                Tp=[]
                for t in T:
                    if not (leq[t,x] or leq[x,t]):
                        Tp.append(t)
                if Tp!=[]:
                    antichains_rec(Ap,Tp[0],Tp[1:],S)
        antichains_rec([],0,range(1,self.size()),S)
        return S

    def subposet_hasse_diagram(self, elements):
        if not isinstance(elements,list):
            raise ValueError, "not a list."
        for element in elements:
            if element not in self:
                raise ValueError, "element not in self"
        relations = []
        for u in elements:
            for v in elements:
                if self.is_lequal(u,v):
                    relations.append([u,v])
        return Poset([elements, relations],cover_relations=False)

