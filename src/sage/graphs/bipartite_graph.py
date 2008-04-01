r"""
Bipartite Graphs

This module implements bipartite graphs.

AUTHORS:
    -- Robert L. Miller (2008-01-20): initial version

TESTS:

    sage: B = graphs.CompleteBipartiteGraph(7,9)
    sage: loads(dumps(B)) == B
    True

"""

#*****************************************************************************
#         Copyright (C) 2008 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

from graph import Graph

class BipartiteGraph(Graph):
    r"""
    Bipartite graph.

    INPUT:
        data -- can be any of the following:
            1.  Empty or None (creates an empty graph).
            2.  An arbitrary graph (finds a bipartition).
            3.  A graph and a bipartition.
            4.  A reduced adjacency matrix.
            5.  A file in alist format.
            6.  From a Networkx bipartite graph.

    A reduced adjacency matrix contains only the non-redundant portion of the
    full adjacency matrix for the bipartite graph.  Specifically, for zero
    matrices of the appropriate size, for the reduced adjacency matrix H, the
    full adjacency matrix is [[0, H'], [H, 0]].

    The alist file format is described at
        http://www.inference.phy.cam.ac.uk/mackay/codes/alist.html

    EXAMPLES:

    1.  No inputs or None for the input creates an empty graph.
        sage: B = BipartiteGraph()
        sage: type(B)
        <class 'sage.graphs.bipartite_graph.BipartiteGraph'>
        sage: B.order()
        0
        sage: B == BipartiteGraph(None)
        True

    2. From a graph: without any more information, finds a bipartition.
        sage: B = BipartiteGraph( graphs.CycleGraph(4) )
        sage: B = BipartiteGraph( graphs.CycleGraph(5) )
        Traceback (most recent call last):
        ...
        TypeError: Input graph is not bipartite!
        sage: G = Graph({0:[5,6], 1:[4,5], 2:[4,6], 3:[4,5,6]})
        sage: B = BipartiteGraph(G)
        sage: B == G
        True

    3. From a graph and a partition. Note that if the input graph is not
    bipartite, then Sage will raise an error. However, if one specifies check =
    False, the offending edges are simply deleted (along with those vertices
    not appearing in either list).  We also lump creating one bipartite graph
    from another into this category.
        sage: P = graphs.PetersenGraph()
        sage: partition = [range(5), range(5,10)]
        sage: B = BipartiteGraph(P, partition)
        Traceback (most recent call last):
        ...
        TypeError: Input graph is not bipartite with respect to the given partition!

        sage: B = BipartiteGraph(P, partition, check=False)
        sage: B.left
        [0, 1, 2, 3, 4]
        sage: B.show()

        sage: G = Graph({0:[5,6], 1:[4,5], 2:[4,6], 3:[4,5,6]})
        sage: B = BipartiteGraph(G)
        sage: B2 = BipartiteGraph(B)
        sage: B == B2
        True
        sage: B3 = BipartiteGraph(G, range(4), range(4,7))
        sage: B3
        Bipartite graph on 7 vertices
        sage: B3 == B2
        True

        sage: G = Graph({0:[], 1:[], 2:[]})
        sage: part = (range(2), [2])
        sage: B = BipartiteGraph(G, part)
        sage: B2 = BipartiteGraph(B)
        sage: B == B2
        True

    4. From a reduced adjacency matrix.
        sage: M = Matrix([(1,1,1,0,0,0,0), (1,0,0,1,1,0,0), \
                          (0,1,0,1,0,1,0), (1,1,0,1,0,0,1)])
        sage: M
        [1 1 1 0 0 0 0]
        [1 0 0 1 1 0 0]
        [0 1 0 1 0 1 0]
        [1 1 0 1 0 0 1]
        sage: H = BipartiteGraph(M); H
        Bipartite graph on 11 vertices
        sage: H.edges()
        [(0, 7, None),
         (0, 8, None),
         (0, 10, None),
         (1, 7, None),
         (1, 9, None),
         (1, 10, None),
         (2, 7, None),
         (3, 8, None),
         (3, 9, None),
         (3, 10, None),
         (4, 8, None),
         (5, 9, None),
         (6, 10, None)]

        sage: M = Matrix([(1, 1, 2, 0, 0), (0, 2, 1, 1, 1), (0, 1, 2, 1, 1)])
        sage: B = BipartiteGraph(M, multiedges=True)
        sage: B.edges()
        [(0, 5, None),
         (1, 5, None),
         (1, 6, None),
         (1, 6, None),
         (1, 7, None),
         (2, 5, None),
         (2, 5, None),
         (2, 6, None),
         (2, 7, None),
         (2, 7, None),
         (3, 6, None),
         (3, 7, None),
         (4, 6, None),
         (4, 7, None)]

         sage: F.<a> = GF(4)
         sage: MS = MatrixSpace(F, 2, 3)
         sage: M = MS.matrix([[0, 1, a+1], [a, 1, 1]])
         sage: B = BipartiteGraph(M, weighted=True)
         sage: B.edges()
         [(0, 4, a), (1, 3, 1), (1, 4, 1), (2, 3, a + 1), (2, 4, 1)]
         sage: B.weighted()
         True

    5. From an alist file.
         sage: file_name = SAGE_TMP + 'deleteme.alist.txt'
         sage: fi = open(file_name, 'w')
         sage: fi.write("7 4 \n 3 4 \n 3 3 1 3 1 1 1 \n 3 3 3 4 \n\
                         1 2 4 \n 1 3 4 \n 1 0 0 \n 2 3 4 \n\
                         2 0 0 \n 3 0 0 \n 4 0 0 \n\
                         1 2 3 0 \n 1 4 5 0 \n 2 4 6 0 \n 1 2 4 7 \n")
         sage: fi.close();
         sage: B = BipartiteGraph(file_name)
         sage: B == H
         True

    6. From a NetworkX bipartite graph.
        sage: import networkx
        sage: G = graphs.OctahedralGraph()
        sage: N = networkx.cliques.make_clique_bipartite(G.networkx_graph())
        sage: B = BipartiteGraph(N)

    """

    def __init__(self, *args, **kwds):
        if len(args) == 0:
            Graph.__init__(self)
            self.left = []; self.right = []
            return
        arg1 = args[0]
        args = args[1:]
        from sage.structure.element import is_Matrix
        if isinstance(arg1, BipartiteGraph):
            Graph.__init__(self, arg1, *args, **kwds)
            self.left, self.right = arg1.left, arg1.right
        elif isinstance(arg1, str):
            Graph.__init__(self, *args, **kwds)
            self.load_afile(arg1)
        elif is_Matrix(arg1):
            # sanity check for mutually exclusive keywords
            if kwds.get('multiedges',False) and kwds.get('weighted',False):
                raise TypeError, "Weighted multi-edge bipartite graphs from reduced adjacency matrix not supported."
            Graph.__init__(self, *args, **kwds)
            ncols = arg1.ncols()
            nrows = arg1.nrows()
            self.left, self.right = range(ncols), range(ncols, nrows+ncols)
            if kwds.get('multiedges',False):
                for ii in range(ncols):
                    for jj in range(nrows):
                        if arg1[jj][ii] != 0:
                            self.add_edges([(ii,jj+ncols)] * arg1[jj][ii])
            elif kwds.get('weighted',False):
                for ii in range(ncols):
                    for jj in range(nrows):
                        if arg1[jj][ii] != 0:
                            self.add_edge((ii, jj+ncols, arg1[jj][ii]))
            else:
                for ii in range(ncols):
                    for jj in range(nrows):
                        if arg1[jj][ii] != 0:
                            self.add_edge((ii, jj+ncols))
        elif isinstance(arg1, Graph) and \
             len(args) > 0 and isinstance(args[0], (list,tuple)) and \
             len(args[0]) == 2 and isinstance(args[0][0], (list,tuple)):
                # Assume that args[0] is a bipartition
                from copy import copy
                left, right = args[0]; left = copy(left); right = copy(right)
                Graph.__init__(self, arg1.subgraph(list(set(left)|set(right))), *args, **kwds)
                if not kwds.has_key('check') or kwds['check']:
                    while len(left) > 0:
                        a = left.pop(0)
                        if len( set( arg1.neighbors(a) ) & set(left) ) != 0:
                            raise TypeError("Input graph is not bipartite with " + \
                             "respect to the given partition!")
                    while len(right) > 0:
                        a = right.pop(0)
                        if len( set( arg1.neighbors(a) ) & set(right) ) != 0:
                            raise TypeError("Input graph is not bipartite with " + \
                             "respect to the given partition!")
                else:
                    while len(left) > 0:
                        a = left.pop(0)
                        a_nbrs = set( arg1.neighbors(a) ) & set(left)
                        if len( a_nbrs ) != 0:
                            self.delete_edges([(a, b) for b in a_nbrs])
                    while len(right) > 0:
                        a = right.pop(0)
                        a_nbrs = set( arg1.neighbors(a) ) & set(right)
                        if len( a_nbrs ) != 0:
                            self.delete_edges([(a, b) for b in a_nbrs])
                self.left, self.right = copy(args[0][0]), copy(args[0][1])
        elif isinstance(arg1, Graph):
            try:
                Graph.__init__(self, arg1, *args, **kwds)
                self.left, self.right = self.bipartite_sets()
                return
            except:
                raise TypeError("Input graph is not bipartite!")
        else:
            import networkx
            Graph.__init__(self, arg1, *args, **kwds)
            if isinstance(arg1, (networkx.XGraph, networkx.Graph)):
                if hasattr(arg1, 'node_type'):
                    # Assume the graph is bipartite
                    self.left = []
                    self.right = []
                    for v in arg1.nodes_iter():
                        if arg1.node_type[v] == 'Bottom':
                            self.left.append(v)
                        elif arg1.node_type[v] == 'Top':
                            self.right.append(v)
                        else:
                            raise TypeError("NetworkX node_type defies bipartite assumtion (is not 'Top' or 'Bottom')")
                else:
                    try:
                        self.left, self.right = self.bipartite_sets()
                    except:
                        raise TypeError("Input graph is not bipartite!")

    def _repr_(self):
        r"""
        Returns a short string representation of self.

        EXAMPLE:
            sage: B = BipartiteGraph(graphs.CycleGraph(16))
            sage: B
            Bipartite cycle graph: graph on 16 vertices

        """
        s = Graph._repr_(self).lower()
        if 'bipartite' in s:
            return s.capitalize()
        else:
            return 'Bipartite ' + s

    def bipartition(self):
        r"""
        Returns the underlying bipartition of the bipartite graph.

        EXAMPLE:
            sage: B = BipartiteGraph( graphs.CycleGraph(4) )
            sage: B.bipartition()
            ([0, 2], [1, 3])

        """
        return (self.left, self.right)

    def project_left(self):
        r"""
        Projects self onto left vertices: edges are 2-paths in the original.

        EXAMPLE:

            sage: B = BipartiteGraph(graphs.CycleGraph(20))
            sage: G = B.project_left()
            sage: G.order(), G.size()
            (10, 10)

        """
        G = Graph()
        G.add_vertices(self.left)
        for v in G:
            for u in self.neighbor_iterator(v):
                G.add_edges([(v,w) for w in self.neighbor_iterator(u)])
        return G

    def project_right(self):
        r"""
        Projects self onto right vertices: edges are 2-paths in the original.

        EXAMPLE:

            sage: B = BipartiteGraph(graphs.CycleGraph(20))
            sage: G = B.project_right()
            sage: G.order(), G.size()
            (10, 10)

        """
        G = Graph()
        G.add_vertices(self.left)
        for v in G:
            for u in self.neighbor_iterator(v):
                G.add_edges([(v,w) for w in self.neighbor_iterator(u)])
        return G

    def plot(self, *args, **kwds):
        r"""
        Overrides Graph's plot function, to illustrate the bipartite nature.

        EXAMPLE:

            sage: B = BipartiteGraph(graphs.CycleGraph(20))
            sage: B.plot()

        """
        if 'pos' not in kwds.keys():
            kwds['pos'] = None
        if kwds['pos'] is None:
            pos = {}
            l_len = len(self.left)
            r_len = len(self.right)
            if l_len == 1:
                pos[self.left[0]] = [-1, 0]
            elif l_len > 1:
                i = 0
                d = 2./(l_len-1)
                for v in self.left:
                    pos[v] = [-1, 1-i*d]
                    i += 1
            if r_len == 1:
                pos[self.right[0]] = [1, 0]
            elif r_len > 1:
                i = 0
                d = 2./(r_len-1)
                for v in self.right:
                    pos[v] = [1, 1-i*d]
                    i += 1
            kwds['pos'] = pos
        return Graph.plot(self, *args, **kwds)

    def load_afile(self, fname):
        r"""
        Loads into the current object the bipartite graph specified in the
        given file name.  This file should follow David MacKay's alist format,
        see
            http://www.inference.phy.cam.ac.uk/mackay/codes/data.html
        for examples and definition of the format.

        EXAMPLE:
            sage: file_name = SAGE_TMP + 'deleteme.alist.txt'
            sage: fi = open(file_name, 'w')
            sage: fi.write("7 4 \n 3 4 \n 3 3 1 3 1 1 1 \n 3 3 3 4 \n\
                            1 2 4 \n 1 3 4 \n 1 0 0 \n 2 3 4 \n\
                            2 0 0 \n 3 0 0 \n 4 0 0 \n\
                            1 2 3 0 \n 1 4 5 0 \n 2 4 6 0 \n 1 2 4 7 \n")
            sage: fi.close();
            sage: B = BipartiteGraph()
            sage: B.load_afile(file_name)
            Bipartite graph on 11 vertices
            sage: B.edges()
            [(0, 7, None),
             (0, 8, None),
             (0, 10, None),
             (1, 7, None),
             (1, 9, None),
             (1, 10, None),
             (2, 7, None),
             (3, 8, None),
             (3, 9, None),
             (3, 10, None),
             (4, 8, None),
             (5, 9, None),
             (6, 10, None)]
             sage: B2 = BipartiteGraph(file_name)
             sage: B2 == B
             True
        """

        # open the file
        try:
            fi = open(fname, 'r')
        except IOError:
            print("Unable to open file <<" + fname + ">>.")
            return None

        # read header information
        num_cols, num_rows = map(int, fi.readline().split())
        max_col_degree, max_row_degree = map(int, fi.readline().split())
        col_degrees = map(int, fi.readline().split())
        row_degrees = map(int, fi.readline().split())

        # sanity checks on header info
        if len(col_degrees) != num_cols:
            print("Invalid Alist format: ")
            print("Number of column degree entries does not match number of columns.")
            return None
        if len(row_degrees) != num_rows:
            print("Invalid Alist format: ")
            print("Number of row degree entries does not match number of rows.")
            return None

        # clear out self
        self.clear()
        self.add_vertices(range(num_cols + num_rows))

        # read adjacency information
        for cidx in range(num_cols):
            for ridx in map(int, fi.readline().split()):
                # A-list uses 1-based indices with 0's as place-holders
                if ridx > 0:
                    self.add_edge(cidx, num_cols + ridx - 1)

        #NOTE:: we could read in the row adjacency information as well to double-check....
        #NOTE:: we could check the actual node degrees against the reported node degrees....

        # now we have all the edges in our graph, just fill in the bipartite partioning
        self.left = list(range(num_cols))
        self.right = list(range(num_cols, num_cols + num_rows))

        # return self for chaining calls if desired
        return self

    def save_afile(self, fname):
        r"""
        Save the graph to file in alist format.

        Saves this graph to file in David MacKay's alist format, see
            http://www.inference.phy.cam.ac.uk/mackay/codes/data.html
        for examples and definition of the format.

        EXAMPLE:
            sage: M = Matrix([(1,1,1,0,0,0,0), (1,0,0,1,1,0,0), \
                              (0,1,0,1,0,1,0), (1,1,0,1,0,0,1)])
            sage: M
            [1 1 1 0 0 0 0]
            [1 0 0 1 1 0 0]
            [0 1 0 1 0 1 0]
            [1 1 0 1 0 0 1]
            sage: b = BipartiteGraph(M)
            sage: file_name = SAGE_TMP + 'deleteme.alist.txt'
            sage: b.save_afile(file_name)
            sage: b2 = BipartiteGraph(file_name)
            sage: b == b2
            True

        TESTS:
            sage: file_name = SAGE_TMP + 'deleteme.alist.txt'
            sage: for order in range(3, 13, 3):
            ...       num_chks = int(order / 3)
            ...       num_vars = order - num_chks
            ...       partition = (range(num_vars), range(num_vars, num_vars+num_chks))
            ...       for idx in range(100):
            ...           g = graphs.RandomGNP(order, 0.5)
            ...           try:
            ...               b = BipartiteGraph(g, partition, check=False)
            ...               b.save_afile(file_name)
            ...               b2 = BipartiteGraph(file_name)
            ...               if b != b2:
            ...                   print "Load/save failed for code with edges:"
            ...                   print b.edges()
            ...           except:
            ...               print "Exception encountered for graph of order "+ str(order)
            ...               print "with edges: "
            ...               g.edges()
            ...               raise

        """

        # open the file
        try:
            fi = open(fname, 'w')
        except IOError:
            print("Unable to open file <<" + fname + ">>.")
            return

        # prep: handy lists, functions for extracting adjacent nodes
        vnodes = self.left
        cnodes = self.right
        max_vdeg = max(self.degree(vnodes))
        max_cdeg = max(self.degree(cnodes))
        num_vnodes = len(vnodes)
        vnbr_str = lambda idx: str(idx - num_vnodes + 1)
        cnbr_str = lambda idx: str(idx + 1)

        # write header information
        fi.write("%d %d\n" % (len(vnodes), len(cnodes)))
        fi.write("%d %d\n" % (max_vdeg, max_cdeg))
        fi.write(" ".join(map(str, self.degree(vnodes))) + "\n")
        fi.write(" ".join(map(str, self.degree(cnodes))) + "\n")
        for vidx in vnodes:
            nbrs = self.neighbors(vidx)
            fi.write(" ".join(map(vnbr_str, nbrs)))
            fi.write(" 0"*(max_vdeg - len(nbrs)) + "\n")
        for cidx in cnodes:
            nbrs = self.neighbors(cidx)
            fi.write(" ".join(map(cnbr_str, nbrs)))
            fi.write(" 0"*(max_cdeg - len(nbrs)) + "\n")

        # done
        fi.close()

        # return self for chaining calls if desired
        return

    def __edge2idx(self, v1, v2):
        r"""
        Translate an edge to its reduced adjacency matrix position.

        Returns (row index, column index) for the given pair of vertices.
        """
        try:
            if v1 in self.left:
                return (self.right.index(v2), self.left.index(v1))
            else:
                return (self.right.index(v1), self.left.index(v2))
        except ValueError:
            raise ValueError("Tried to map invalid edge (%d,%d) to vertex indices" \
                                 % (v1, v2))

    def reduced_adjacency_matrix(self, sparse=True):
        r"""
        Return the reduced adjacency matrix for the given graph.

        A reduced adjacency matrix contains only the non-redundant portion of the
        full adjacency matrix for the bipartite graph.  Specifically, for zero
        matrices of the appropriate size, for the reduced adjacency matrix H, the
        full adjacency matrix is [[0, H'], [H, 0]].

        This method supports the named argument 'sparse' which defaults to
        True.  When enabled, the returned matrix will be sparse.

        EXAMPLES:

        Bipartite graphs that are not weighted will return a matrix over ZZ.
            sage: M = Matrix([(1,1,1,0,0,0,0), (1,0,0,1,1,0,0), \
                              (0,1,0,1,0,1,0), (1,1,0,1,0,0,1)])
            sage: B = BipartiteGraph(M)
            sage: N = B.reduced_adjacency_matrix()
            sage: N
            [1 1 1 0 0 0 0]
            [1 0 0 1 1 0 0]
            [0 1 0 1 0 1 0]
            [1 1 0 1 0 0 1]
            sage: N == M
            True
            sage: N[0,0].parent()
            Integer Ring

        Multi-edge graphs also return a matrix over ZZ.
            sage: M = Matrix([(1, 1, 2, 0, 0), (0, 2, 1, 1, 1), (0, 1, 2, 1, 1)])
            sage: B = BipartiteGraph(M, multiedges=True)
            sage: N = B.reduced_adjacency_matrix()
            sage: N == M
            True
            sage: N[0,0].parent()
            Integer Ring

        Weighted graphs will return a matrix over the ring given by their
        (first) weights.
            sage: F.<a> = GF(4)
            sage: MS = MatrixSpace(F, 2, 3)
            sage: M = MS.matrix([[0, 1, a+1], [a, 1, 1]])
            sage: B = BipartiteGraph(M, weighted=True)
            sage: N = B.reduced_adjacency_matrix(sparse=False)
            sage: N == M
            True
            sage: N[0,0].parent()
            Finite Field in a of size 2^2

        TESTS:
            sage: B = BipartiteGraph()
            sage: B.reduced_adjacency_matrix()
            []
            sage: M = Matrix([[0,0], [0,0]])
            sage: BipartiteGraph(M).reduced_adjacency_matrix() == M
            True
            sage: M = Matrix([[10,2/3], [0,0]])
            sage: B = BipartiteGraph(M, weighted=True)
            sage: M == B.reduced_adjacency_matrix()
            True
        """
        if self.multiple_edges() and self.weighted():
            raise NotImplementedError, "Don't know how to represent weights for a multigraph."
        if self.is_directed():
            raise NotImplementedError, "Reduced adjacency matrix does not exist for directed graphs."

        # create dictionary of edges, values are weights for weighted graph,
        # otherwise the number of edges (always 1 for simple graphs)
        D = {}
        if self.weighted():
            for (v1, v2, weight) in self.edge_iterator():
                D[self.__edge2idx(v1,v2)] = weight
        else:
            # if we're normal or multi-edge, just create the matrix over ZZ
            for (v1, v2, name) in self.edge_iterator():
                idx = self.__edge2idx(v1, v2)
                if D.has_key(idx):
                    D[idx] = 1 + D[idx]
                else:
                    D[idx] = 1

        # now construct and return the matrix from the dictionary we created
        from sage.matrix.constructor import matrix
        return matrix(len(self.right), len(self.left), D, sparse=sparse)
