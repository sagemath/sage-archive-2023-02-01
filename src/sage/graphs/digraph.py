r"""
Directed graphs

This module implements functions and operations involving directed
graphs.

Class and methods
-----------------
"""

from sage.rings.integer import Integer
from sage.misc.misc import deprecated_function_alias
import sage.graphs.generic_graph_pyx as generic_graph_pyx
from sage.graphs.generic_graph import GenericGraph

class DiGraph(GenericGraph):
    """
    Directed graph.

    A digraph or directed graph is a set of vertices connected by oriented
    edges (cf. http://en.wikipedia.org/wiki/Digraph_%28mathematics%29 ).

    One can very easily create a directed graph in Sage by typing::

        sage: g = DiGraph()

    By typing the name of the digraph, one can get some basic information
    about it::

        sage: g
        Digraph on 0 vertices

    This digraph is not very interesting as it is by default the empty
    graph. But Sage contains several pre-defined digraph classes that can
    be listed this way:

    * Within a Sage sessions, type ``digraphs.``
      (Do not press "Enter", and do not forget the final period "." )
    * Hit "tab".

    You will see a list of methods which will construct named digraphs. For
    example::

        sage: g = digraphs.ButterflyGraph(3)
        sage: g.plot()

    You can also use the collection of pre-defined graphs, then create a
    digraph from them. ::

        sage: g = DiGraph(graphs.PetersenGraph())
        sage: g.plot()

    Calling ``Digraph`` on a graph returns the original graph in which every
    edge is replaced by two different edges going toward opposite directions.

    In order to obtain more information about these digraph constructors,
    access the documentation by typing ``digraphs.RandomDirectedGNP?``.

    Once you have defined the digraph you want, you can begin to work on it
    by using the almost 200 functions on graphs and digraphs in the Sage
    library! If your digraph is named ``g``, you can list these functions as
    previously this way

    * Within a Sage session, type ``g.``
      (Do not press "Enter", and do not forget the final period "." )
    * Hit "tab".

    As usual, you can get some information about what these functions do by
    typing (e.g. if you want to know about the ``diameter()`` method)
    ``g.diameter?``.

    If you have defined a digraph ``g`` having several connected components
    ( i.e. ``g.is_connected()`` returns False ), you can print each one of its
    connected components with only two lines::

        sage: for component in g.connected_components():
        ...      g.subgraph(component).plot()

    The same methods works for strongly connected components ::
        sage: for component in g.strongly_connected_components():
        ...      g.subgraph(component).plot()


    INPUT:

    -  ``data`` -  can be any of the following:

       #.  A dictionary of dictionaries

       #.  A dictionary of lists

       #.  A numpy matrix or ndarray

       #.  A Sage adjacency matrix or incidence matrix

       #.  A pygraphviz agraph

       #.  A SciPy sparse matrix

       #.  A NetworkX digraph

    -  ``pos`` - a positioning dictionary: for example, the
       spring layout from NetworkX for the 5-cycle is::

         {0: [-0.91679746, 0.88169588],
          1: [ 0.47294849, 1.125     ],
          2: [ 1.125     ,-0.12867615],
          3: [ 0.12743933,-1.125     ],
          4: [-1.125     ,-0.50118505]}

    -  ``name`` - (must be an explicitly named parameter,
       i.e., name="complete") gives the graph a name

    -  ``loops`` - boolean, whether to allow loops (ignored
       if data is an instance of the DiGraph class)

    -  ``multiedges`` - boolean, whether to allow multiple
       edges (ignored if data is an instance of the DiGraph class)

    -  ``weighted`` - whether digraph thinks of itself as
       weighted or not. See self.weighted()

    -  ``format`` - if None, DiGraph tries to guess- can be
       several values, including:

       -  ``'adjacency_matrix'`` - a square Sage matrix M,
          with M[i,j] equal to the number of edges {i,j}

       -  ``'incidence_matrix'`` - a Sage matrix, with one
          column C for each edge, where if C represents {i, j}, C[i] is -1
          and C[j] is 1

       -  ``'weighted_adjacency_matrix'`` - a square Sage
          matrix M, with M[i,j] equal to the weight of the single edge {i,j}.
          Given this format, weighted is ignored (assumed True).

    -  ``boundary`` - a list of boundary vertices, if none,
       digraph is considered as a 'digraph without boundary'

    -  ``implementation`` - what to use as a backend for
       the graph. Currently, the options are either 'networkx' or
       'c_graph'

    -  ``sparse`` - only for implementation == 'c_graph'.
       Whether to use sparse or dense graphs as backend. Note that
       currently dense graphs do not have edge labels, nor can they be
       multigraphs

    -  ``vertex_labels`` - only for implementation == 'c_graph'.
       Whether to allow any object as a vertex (slower), or
       only the integers 0, ..., n-1, where n is the number of vertices.


    EXAMPLES:


    #. A dictionary of dictionaries::

            sage: g = DiGraph({0:{1:'x',2:'z',3:'a'}, 2:{5:'out'}}); g
            Digraph on 5 vertices

       The labels ('x', 'z', 'a', 'out') are labels for edges. For
       example, 'out' is the label for the edge from 2 to 5. Labels can be
       used as weights, if all the labels share some common parent.

    #. A dictionary of lists::

            sage: g = DiGraph({0:[1,2,3], 2:[4]}); g
            Digraph on 5 vertices

    #. A list of vertices and a function describing adjacencies. Note
       that the list of vertices and the function must be enclosed in a
       list (i.e., [list of vertices, function]).

       We construct a graph on the integers 1 through 12 such that there
       is a directed edge from i to j if and only if i divides j.

       ::

            sage: g=DiGraph([[1..12],lambda i,j: i!=j and i.divides(j)])
            sage: g.vertices()
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
            sage: g.adjacency_matrix()
            [0 1 1 1 1 1 1 1 1 1 1 1]
            [0 0 0 1 0 1 0 1 0 1 0 1]
            [0 0 0 0 0 1 0 0 1 0 0 1]
            [0 0 0 0 0 0 0 1 0 0 0 1]
            [0 0 0 0 0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 1]
            [0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0]

    #. A numpy matrix or ndarray::

            sage: import numpy
            sage: A = numpy.array([[0,1,0],[1,0,0],[1,1,0]])
            sage: DiGraph(A)
            Digraph on 3 vertices

    #. A Sage matrix: Note: If format is not specified, then Sage
       assumes a square matrix is an adjacency matrix, and a nonsquare
       matrix is an incidence matrix.

       - an adjacency matrix::

            sage: M = Matrix([[0, 1, 1, 1, 0],[0, 0, 0, 0, 0],[0, 0, 0, 0, 1],[0, 0, 0, 0, 0],[0, 0, 0, 0, 0]]); M
            [0 1 1 1 0]
            [0 0 0 0 0]
            [0 0 0 0 1]
            [0 0 0 0 0]
            [0 0 0 0 0]
            sage: DiGraph(M)
            Digraph on 5 vertices

       - an incidence matrix::

            sage: M = Matrix(6, [-1,0,0,0,1, 1,-1,0,0,0, 0,1,-1,0,0, 0,0,1,-1,0, 0,0,0,1,-1, 0,0,0,0,0]); M
            [-1  0  0  0  1]
            [ 1 -1  0  0  0]
            [ 0  1 -1  0  0]
            [ 0  0  1 -1  0]
            [ 0  0  0  1 -1]
            [ 0  0  0  0  0]
            sage: DiGraph(M)
            Digraph on 6 vertices

    #. A c_graph implemented DiGraph can be constructed from a
       networkx implemented DiGraph::

            sage: D = DiGraph({0:[1],1:[2],2:[0]}, implementation="networkx")
            sage: E = DiGraph(D,implementation="c_graph")
            sage: D == E
            True

    #. A dig6 string: Sage automatically recognizes whether a string is
       in dig6 format, which is a directed version of graph6::

            sage: D = DiGraph('IRAaDCIIOWEOKcPWAo')
            sage: D
            Digraph on 10 vertices

            sage: D = DiGraph('IRAaDCIIOEOKcPWAo')
            Traceback (most recent call last):
            ...
            RuntimeError: The string (IRAaDCIIOEOKcPWAo) seems corrupt: for n = 10, the string is too short.

            sage: D = DiGraph("IRAaDCI'OWEOKcPWAo")
            Traceback (most recent call last):
            ...
            RuntimeError: The string seems corrupt: valid characters are
            ?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~

    #. A NetworkX XDiGraph::

            sage: import networkx
            sage: g = networkx.MultiDiGraph({0:[1,2,3], 2:[4]})
            sage: DiGraph(g)
            Digraph on 5 vertices


    #. A NetworkX digraph::

            sage: import networkx
            sage: g = networkx.DiGraph({0:[1,2,3], 2:[4]})
            sage: DiGraph(g)
            Digraph on 5 vertices

       Note that in all cases, we copy the NetworkX structure.

       ::

            sage: import networkx
            sage: g = networkx.DiGraph({0:[1,2,3], 2:[4]})
            sage: G = DiGraph(g, implementation='networkx')
            sage: H = DiGraph(g, implementation='networkx')
            sage: G._backend._nxg is H._backend._nxg
            False


    """
    _directed = True

    def __init__(self, data=None, pos=None, loops=None, format=None,
                 boundary=[], weighted=None, implementation='c_graph',
                 sparse=True, vertex_labels=True, **kwds):
        """
        TESTS::

            sage: D = DiGraph()
            sage: loads(dumps(D)) == D
            True

            sage: a = matrix(2,2,[1,2,0,1])
            sage: DiGraph(a,sparse=True).adjacency_matrix() == a
            True

            sage: a = matrix(2,2,[3,2,0,1])
            sage: DiGraph(a,sparse=True).adjacency_matrix() == a
            True

        The positions are copied when the DiGraph is built from
        another DiGraph or from a Graph ::

            sage: g = DiGraph(graphs.PetersenGraph())
            sage: h = DiGraph(g)
            sage: g.get_pos() == h.get_pos()
            True
            sage: g.get_pos() == graphs.PetersenGraph().get_pos()
            True
        """
        from sage.all import walltime
        msg = ''
        GenericGraph.__init__(self)
        multiedges = kwds.get('multiedges', None)
        from sage.structure.element import is_Matrix
        from sage.misc.misc import uniq
        if format is None and isinstance(data, str):
            format = 'dig6'
            if data[:8] == ">>dig6<<":
                data = data[8:]
        if format is None and is_Matrix(data):
            if data.is_square():
                format = 'adjacency_matrix'
            else:
                format = 'incidence_matrix'
                msg += "Non-symmetric or non-square matrix assumed to be an incidence matrix: "
        if format is None and isinstance(data, DiGraph):
            format = 'DiGraph'
        from sage.graphs.all import Graph
        if format is None and isinstance(data, Graph):
            data = data.to_directed()
            format = 'DiGraph'
        if format is None and isinstance(data,list) and \
           len(data)>=2 and callable(data[1]):
            format = 'rule'
        if format is None and isinstance(data,dict):
            keys = data.keys()
            if len(keys) == 0: format = 'dict_of_dicts'
            else:
                if isinstance(data[keys[0]], list):
                    format = 'dict_of_lists'
                elif isinstance(data[keys[0]], dict):
                    format = 'dict_of_dicts'
        if format is None and hasattr(data, 'adj'):
            import networkx
            if isinstance(data, (networkx.Graph, networkx.MultiGraph)):
                data = data.to_directed()
                format = 'NX'
            elif isinstance(data, (networkx.DiGraph, networkx.MultiDiGraph)):
                format = 'NX'
        if format is None and isinstance(data, (int, Integer)):
            format = 'int'
        if format is None and data is None:
            format = 'int'
            data = 0
        if format is None:
            import networkx
            data = networkx.MultiDiGraph(data)
            format = 'NX'

        # At this point, format has been set.
        verts = None

        if format == 'dig6':
            if loops      is None: loops      = False
            if weighted   is None: weighted   = False
            if multiedges is None: multiedges = False
            if not isinstance(data, str):
                raise ValueError('If input format is dig6, then data must be a string.')
            n = data.find('\n')
            if n == -1:
                n = len(data)
            ss = data[:n]
            n, s = generic_graph_pyx.N_inverse(ss)
            m = generic_graph_pyx.D_inverse(s, n)
            expected = n**2
            if len(m) > expected:
                raise RuntimeError("The string (%s) seems corrupt: for n = %d, the string is too long."%(ss,n))
            elif len(m) < expected:
                raise RuntimeError("The string (%s) seems corrupt: for n = %d, the string is too short."%(ss,n))
            num_verts = n
        elif format in ['adjacency_matrix', 'incidence_matrix']:
            assert is_Matrix(data)
        if format == 'weighted_adjacency_matrix':
            if weighted is False:
                raise ValueError("Format was weighted_adjacency_matrix but weighted was False.")
            if weighted   is None: weighted   = True
            if multiedges is None: multiedges = False
            format = 'adjacency_matrix'
        if format == 'adjacency_matrix':
            entries = uniq(data.list())
            for e in entries:
                try:
                    e = int(e)
                    assert e >= 0
                except:
                    if weighted is False:
                        raise ValueError("Non-weighted digraph's"+
                        " adjacency matrix must have only nonnegative"+
                        " integer entries")
                    weighted = True
                    if multiedges is None: multiedges = False
                    break
            if multiedges is None: multiedges = (sorted(entries) != [0,1])
            if weighted is None: weighted = False
            for i in xrange(data.nrows()):
                if data[i,i] != 0:
                    if loops is None: loops = True
                    elif not loops:
                        raise ValueError("Non-looped digraph's adjacency"+
                        " matrix must have zeroes on the diagonal.")
                    break
            num_verts = data.nrows()
        elif format == 'incidence_matrix':
            try:
                positions = []
                for c in data.columns():
                    NZ = c.nonzero_positions()
                    if len(NZ) != 2:
                        msg += "There must be two nonzero entries (-1 & 1) per column."
                        assert False
                    L = uniq(c.list())
                    L.sort()
                    if L != [-1,0,1]:
                        msg += "Each column represents an edge: -1 goes to 1."
                        assert False
                    if c[NZ[0]] == -1:
                        positions.append(tuple(NZ))
                    else:
                        positions.append((NZ[1],NZ[0]))
                if loops      is None: loops     = False
                if weighted   is None: weighted  = False
                if multiedges is None:
                    total = len(positions)
                    multiedges = (  len(uniq(positions)) < total  )
            except AssertionError:
                raise ValueError(msg)
            num_verts = data.nrows()
        elif format == 'DiGraph':
            if loops is None: loops = data.allows_loops()
            elif not loops and data.has_loops():
                raise ValueError("No loops but input digraph has loops.")
            if multiedges is None: multiedges = data.allows_multiple_edges()
            elif not multiedges:
                e = data.edges(labels=False)
                if len(e) != len(uniq(e)):
                    raise ValueError("No multiple edges but input digraph"+
                    " has multiple edges.")
            if weighted is None: weighted = data.weighted()
            num_verts = data.num_verts()
            verts = data.vertex_iterator()
            if data.get_pos() is not None:
                pos = data.get_pos().copy()
        elif format == 'rule':
            f = data[1]
            if loops is None: loops = any(f(v,v) for v in data[0])
            if multiedges is None: multiedges = False
            if weighted is None: weighted = False
            num_verts = len(data[0])
            verts = data[0]
        elif format == 'dict_of_dicts':
            if not all(isinstance(data[u], dict) for u in data):
                raise ValueError("Input dict must be a consistent format.")
            verts = set(data.keys())
            if loops is None or loops is False:
                for u in data:
                    if u in data[u]:
                        if loops is None: loops = True
                        elif loops is False:
                            raise ValueError("No loops but dict has loops.")
                        break
                if loops is None: loops = False
            if weighted is None: weighted = False
            for u in data:
                for v in data[u]:
                    if v not in verts: verts.add(v)
                    if multiedges is not False and not isinstance(data[u][v], list):
                        if multiedges is None: multiedges = False
                        if multiedges:
                            raise ValueError("Dict of dicts for multidigraph must be in the format {v : {u : list}}")
            if multiedges is None: multiedges = True
            num_verts = len(verts)
        elif format == 'dict_of_lists':
            if not all(isinstance(data[u], list) for u in data):
                raise ValueError("Input dict must be a consistent format.")
            verts = set(data.keys())
            if loops is None or loops is False:
                for u in data:
                    if u in data[u]:
                        if loops is None: loops = True
                        elif loops is False:
                            raise ValueError("No loops but dict has loops.")
                        break
                if loops is None: loops = False
            if weighted is None: weighted = False
            for u in data:
                verts = verts.union([v for v in data[u] if v not in verts])
                if len(uniq(data[u])) != len(data[u]):
                    if multiedges is False:
                        raise ValueError("Non-multidigraph input dict has multiple edges (%s,%s)"%(u, choice([v for v in data[u] if data[u].count(v) > 1])))
                    if multiedges is None: multiedges = True
            if multiedges is None: multiedges = False
            num_verts = len(verts)
        elif format == 'NX':
            if weighted is None:
                if isinstance(data, networkx.DiGraph):
                    weighted = False
                    if multiedges is None:
                        multiedges = False
                    if loops is None:
                        loops = False
                else:
                    weighted = True
                    if multiedges is None:
                        multiedges = data.multiedges
                    if loops is None:
                        loops = data.selfloops
            num_verts = data.order()
            verts = data.nodes()
            data = data.adj
            format = 'dict_of_dicts'
        elif format == 'int':
            if weighted   is None: weighted   = False
            if multiedges is None: multiedges = False
            if loops      is None: loops      = False
            num_verts = data

        # weighted, multiedges, loops, verts and num_verts should now be set

        if implementation == 'networkx':
            import networkx
            from sage.graphs.base.graph_backends import NetworkXGraphBackend
            if format == 'DiGraph':
                self._backend = NetworkXGraphBackend(data.networkx_graph())
                self._weighted = weighted
                self.allow_loops(loops)
                self.allow_multiple_edges(multiedges)
            else:
                self._backend = NetworkXGraphBackend(networkx.MultiDiGraph())
                self._weighted = weighted
                self.allow_loops(loops)
                self.allow_multiple_edges(multiedges)
                if verts is not None:
                    self.add_vertices(verts)
                else:
                    self.add_vertices(range(num_verts))
        elif implementation == 'c_graph':
            if multiedges or weighted:
                if not sparse:
                    raise RuntimeError("Multiedge and weighted c_graphs must be sparse.")
            from sage.graphs.base.sparse_graph import SparseGraphBackend
            from sage.graphs.base.dense_graph import DenseGraphBackend
            CGB = SparseGraphBackend if sparse else DenseGraphBackend
            if format == 'DiGraph':
                self._backend = CGB(0, directed=True)
                self.add_vertices(verts)
                self._weighted = weighted
                self.allow_loops(loops, check=False)
                self.allow_multiple_edges(multiedges, check=False)
                for u,v,l in data.edge_iterator():
                    self._backend.add_edge(u,v,l,True)
            else:
                self._backend = CGB(0, directed=True)
                if verts is not None:
                    self._backend = CGB(0, directed=True)
                    self.add_vertices(verts)
                else:
                    self._backend = CGB(num_verts, directed=True)
                self._weighted = weighted
                self.allow_loops(loops, check=False)
                self.allow_multiple_edges(multiedges, check=False)
            self._backend.directed = True
        else:
            raise NotImplementedError("Supported implementations: networkx, c_graph.")

        if format == 'dig6':
            k = 0
            for i in xrange(n):
                for j in xrange(n):
                    if m[k] == '1':
                        self.add_edge(i, j)
                    k += 1
        elif format == 'adjacency_matrix':
            e = []
            if weighted:
                for i,j in data.nonzero_positions():
                    e.append((i,j,data[i][j]))
            elif multiedges:
                for i,j in data.nonzero_positions():
                    e += [(i,j)]*int(data[i][j])
            else:
                for i,j in data.nonzero_positions():
                    e.append((i,j))
            self.add_edges(e)
        elif format == 'incidence_matrix':
            self.add_edges(positions)
        elif format == 'DiGraph':
            self.name(data.name())
        elif format == 'rule':
            verts = list(verts)
            for u in xrange(num_verts):
                for v in xrange(num_verts):
                    uu,vv = verts[u], verts[v]
                    if f(uu,vv):
                        self.add_edge(uu,vv)
        elif format == 'dict_of_dicts':
            for u in data:
                for v in data[u]:
                    if multiedges:
                        self.add_edges([(u,v,l) for l in data[u][v]])
                    else:
                        self.add_edge((u,v,data[u][v]))
        elif format == 'dict_of_lists':
            for u in data:
                for v in data[u]:
                    self.add_edge(u,v)
        else:
            assert format == 'int'
        self._pos = pos
        self._boundary = boundary
        name = kwds.get('name', None)
        if format != 'DiGraph' or name is not None:
            self.name(name)

    ### Formats

    def dig6_string(self):
        """
        Returns the dig6 representation of the digraph as an ASCII string.
        Valid for single (no multiple edges) digraphs on 0 to 262143
        vertices.

        EXAMPLES::

            sage: D = DiGraph()
            sage: D.dig6_string()
            '?'
            sage: D.add_edge(0,1)
            sage: D.dig6_string()
            'AO'
        """
        n = self.order()
        if n > 262143:
            raise ValueError('dig6 format supports graphs on 0 to 262143 vertices only.')
        elif self.has_multiple_edges():
            raise ValueError('dig6 format does not support multiple edges.')
        else:
            return generic_graph_pyx.N(n) + generic_graph_pyx.R(self._bit_vector())

    ### Attributes

    def is_directed(self):
        """
        Since digraph is directed, returns True.

        EXAMPLES::

            sage: DiGraph().is_directed()
            True
        """
        return True

    ### Properties

    def is_directed_acyclic(self):
        """
        Returns whether the digraph is acyclic or not.

        A directed graph is acyclic if for any vertex v, there is no
        directed path that starts and ends at v. Every directed acyclic
        graph (dag) corresponds to a partial ordering of its vertices,
        however multiple dags may lead to the same partial ordering.

        EXAMPLES::

            sage: D = DiGraph({ 0:[1,2,3], 4:[2,5], 1:[8], 2:[7], 3:[7], 5:[6,7], 7:[8], 6:[9], 8:[10], 9:[10] })
            sage: D.plot(layout='circular').show()
            sage: D.is_directed_acyclic()
            True

        ::

            sage: D.add_edge(9,7)
            sage: D.is_directed_acyclic()
            True

        ::

            sage: D.add_edge(7,4)
            sage: D.is_directed_acyclic()
            False
        """
        try:
            S = self.topological_sort()
            return True
        except:
            return False

    def to_directed(self):
        """
        Since the graph is already directed, simply returns a copy of
        itself.

        EXAMPLES::

            sage: DiGraph({0:[1,2,3],4:[5,1]}).to_directed()
            Digraph on 6 vertices
        """
        from copy import copy
        return copy(self)

    def to_undirected(self, implementation='c_graph', sparse=None):
        """
        Returns an undirected version of the graph. Every directed edge
        becomes an edge.

        EXAMPLES::

            sage: D = DiGraph({0:[1,2],1:[0]})
            sage: G = D.to_undirected()
            sage: D.edges(labels=False)
            [(0, 1), (0, 2), (1, 0)]
            sage: G.edges(labels=False)
            [(0, 1), (0, 2)]
        """
        if sparse is None:
            from sage.graphs.base.dense_graph import DenseGraphBackend
            sparse = (not isinstance(self._backend, DenseGraphBackend))
        from sage.graphs.all import Graph
        G = Graph(name=self.name(), pos=self._pos, boundary=self._boundary,
                  multiedges=self.allows_multiple_edges(), loops=self.allows_loops(),
                  implementation=implementation, sparse=sparse)
        G.name(self.name())
        G.add_vertices(self.vertex_iterator())
        G.add_edges(self.edge_iterator())
        if hasattr(self, '_embedding'):
            from copy import copy
            G._embedding = copy(self._embedding)
        G._weighted = self._weighted
        return G

    ### Edge Handlers

    def incoming_edge_iterator(self, vertices=None, labels=True):
        """
        Return an iterator over all arriving edges from vertices, or over
        all edges if vertices is None.

        EXAMPLES::

            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: for a in D.incoming_edge_iterator([0]):
            ...    print a
            (1, 0, None)
            (4, 0, None)
        """
        if vertices is None:
            vertices = self
        elif vertices in self:
            vertices = [vertices]
        else:
            vertices = [v for v in vertices if v in self]
        return self._backend.iterator_in_edges(vertices, labels)

    def incoming_edges(self, vertices=None, labels=True):
        """
        Returns a list of edges arriving at vertices.

        INPUT:


        -  ``labels`` - if False, each edge is a tuple (u,v) of
           vertices.


        EXAMPLES::

            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.incoming_edges([0])
            [(1, 0, None), (4, 0, None)]
        """
        return list(self.incoming_edge_iterator(vertices, labels=labels))

    def outgoing_edge_iterator(self, vertices=None, labels=True):
        """
        Return an iterator over all departing edges from vertices, or over
        all edges if vertices is None.

        EXAMPLES::

            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: for a in D.outgoing_edge_iterator([0]):
            ...    print a
            (0, 1, None)
            (0, 2, None)
            (0, 3, None)
        """
        if vertices is None:
            vertices = self
        elif vertices in self:
            vertices = [vertices]
        else:
            vertices = [v for v in vertices if v in self]
        return self._backend.iterator_out_edges(vertices, labels)

    def outgoing_edges(self, vertices=None, labels=True):
        """
        Returns a list of edges departing from vertices.

        INPUT:


        -  ``labels`` - if False, each edge is a tuple (u,v) of
           vertices.


        EXAMPLES::

            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.outgoing_edges([0])
            [(0, 1, None), (0, 2, None), (0, 3, None)]
        """
        return list(self.outgoing_edge_iterator(vertices, labels=labels))

    def neighbor_in_iterator(self, vertex):
        """
        Returns an iterator over the in-neighbors of vertex.

        An vertex `u` is an in-neighbor of a vertex `v` if `uv` in an edge.

        EXAMPLES::

            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: for a in D.neighbor_in_iterator(0):
            ...    print a
            1
            4
        """
        return iter(set(self._backend.iterator_in_nbrs(vertex)))

    predecessor_iterator = deprecated_function_alias(neighbor_in_iterator, 'Sage Version 4.3')

    def neighbors_in(self, vertex):
        """
        Returns the list of the in-neighbors of a given vertex.

        An vertex `u` is an in-neighbor of a vertex `v` if `uv` in an edge.

        EXAMPLES::

            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.neighbors_in(0)
            [1, 4]
        """
        return list(self.neighbor_in_iterator(vertex))

    predecessors = deprecated_function_alias(neighbors_in, 'Sage Version 4.3')

    def neighbor_out_iterator(self, vertex):
        """
        Returns an iterator over the out-neighbors of a given vertex.

        An vertex `u` is an out-neighbor of a vertex `v` if `vu` in an edge.

        EXAMPLES::

            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: for a in D.neighbor_out_iterator(0):
            ...    print a
            1
            2
            3
        """
        return iter(set(self._backend.iterator_out_nbrs(vertex)))

    successor_iterator = deprecated_function_alias(neighbor_out_iterator, 'Sage Version 4.3')

    def neighbors_out(self, vertex):
        """
        Returns the list of the out-neighbors of a given vertex.

        An vertex `u` is an out-neighbor of a vertex `v` if `vu` in an edge.

        EXAMPLES::

            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.neighbors_out(0)
            [1, 2, 3]
        """
        return list(self.neighbor_out_iterator(vertex))

    successors = deprecated_function_alias(neighbors_out, 'Sage Version 4.3')

    ### Degree functions

    def in_degree(self, vertices=None, labels=False):
        """
        Same as degree, but for in degree.

        EXAMPLES::

            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.in_degree(vertices = [0,1,2], labels=True)
            {0: 2, 1: 2, 2: 2}
            sage: D.in_degree()
            [2, 2, 2, 2, 1, 1]
            sage: G = graphs.PetersenGraph().to_directed()
            sage: G.in_degree(0)
            3
        """
        if vertices in self:
            for d in self.in_degree_iterator(vertices):
                return d # (weird, but works: only happens once!)
        elif labels:
            di = {}
            for v, d in self.in_degree_iterator(vertices, labels=labels):
                di[v] = d
            return di
        else:
            return list(self.in_degree_iterator(vertices, labels=labels))

    def in_degree_iterator(self, vertices=None, labels=False):
        """
        Same as degree_iterator, but for in degree.

        EXAMPLES::

            sage: D = graphs.Grid2dGraph(2,4).to_directed()
            sage: for i in D.in_degree_iterator():
            ...    print i
            3
            3
            2
            3
            2
            2
            2
            3
            sage: for i in D.in_degree_iterator(labels=True):
            ...    print i
            ((0, 1), 3)
            ((1, 2), 3)
            ((0, 0), 2)
            ((0, 2), 3)
            ((1, 3), 2)
            ((1, 0), 2)
            ((0, 3), 2)
            ((1, 1), 3)
        """
        if vertices is None:
            vertices = self.vertex_iterator()
        elif vertices in self:
            d = 0
            for e in self.incoming_edge_iterator(vertices):
                d += 1
            yield d
            return
        if labels:
            for v in vertices:
                d = 0
                for e in self.incoming_edge_iterator(v):
                    d += 1
                yield (v, d)
        else:
            for v in vertices:
                d = 0
                for e in self.incoming_edge_iterator(v):
                    d += 1
                yield d

    def in_degree_sequence(self):
        r"""
        Return the indegree sequence of this digraph.

        EXAMPLES:

        The indegree sequences of two digraphs::

            sage: g = DiGraph({1: [2, 5, 6], 2: [3, 6], 3: [4, 6], 4: [6], 5: [4, 6]})
            sage: g.in_degree_sequence()
            [5, 2, 1, 1, 1, 0]

        ::

            sage: V = [2, 3, 5, 7, 8, 9, 10, 11]
            sage: E = [[], [8, 10], [11], [8, 11], [9], [], [], [2, 9, 10]]
            sage: g = DiGraph(dict(zip(V, E)))
            sage: g.in_degree_sequence()
            [2, 2, 2, 2, 1, 0, 0, 0]
        """
        return sorted(self.in_degree_iterator(), reverse=True)

    def out_degree(self, vertices=None, labels=False):
        """
        Same as degree, but for out degree.

        EXAMPLES::

            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.out_degree(vertices = [0,1,2], labels=True)
            {0: 3, 1: 2, 2: 1}
            sage: D.out_degree()
            [3, 2, 1, 1, 2, 1]
        """
        if vertices in self:
            for d in self.out_degree_iterator(vertices):
                return d # (weird, but works: only happens once!)
        elif labels:
            di = {}
            for v, d in self.out_degree_iterator(vertices, labels=labels):
                di[v] = d
            return di
        else:
            return list(self.out_degree_iterator(vertices, labels=labels))

    def out_degree_iterator(self, vertices=None, labels=False):
        """
        Same as degree_iterator, but for out degree.

        EXAMPLES::

            sage: D = graphs.Grid2dGraph(2,4).to_directed()
            sage: for i in D.out_degree_iterator():
            ...    print i
            3
            3
            2
            3
            2
            2
            2
            3
            sage: for i in D.out_degree_iterator(labels=True):
            ...    print i
            ((0, 1), 3)
            ((1, 2), 3)
            ((0, 0), 2)
            ((0, 2), 3)
            ((1, 3), 2)
            ((1, 0), 2)
            ((0, 3), 2)
            ((1, 1), 3)
        """
        if vertices is None:
            vertices = self.vertex_iterator()
        elif vertices in self:
            d = 0
            for e in self.outgoing_edge_iterator(vertices):
                d += 1
            yield d
            return
        if labels:
            for v in vertices:
                d = 0
                for e in self.outgoing_edge_iterator(v):
                    d += 1
                yield (v, d)
        else:
            for v in vertices:
                d = 0
                for e in self.outgoing_edge_iterator(v):
                    d += 1
                yield d

    def out_degree_sequence(self):
        r"""
        Return the outdegree sequence of this digraph.

        EXAMPLES:

        The outdegree sequences of two digraphs::

            sage: g = DiGraph({1: [2, 5, 6], 2: [3, 6], 3: [4, 6], 4: [6], 5: [4, 6]})
            sage: g.out_degree_sequence()
            [3, 2, 2, 2, 1, 0]

        ::

            sage: V = [2, 3, 5, 7, 8, 9, 10, 11]
            sage: E = [[], [8, 10], [11], [8, 11], [9], [], [], [2, 9, 10]]
            sage: g = DiGraph(dict(zip(V, E)))
            sage: g.out_degree_sequence()
            [3, 2, 2, 1, 1, 0, 0, 0]
        """
        return sorted(self.out_degree_iterator(), reverse=True)


    def feedback_edge_set(self,value_only=False):
        r"""
        Computes the minimum feedback edge set of a digraph
        ( also called feedback arc set ).

        The minimum feedback edge set of a digraph is a set of edges
        that intersect all the circuits of the digraph.
        Equivalently, a minimum feedback arc set of a DiGraph is a set
        `S` of arcs such that the digraph `G-S` is acyclic.

        For more informations, see
        ( http://en.wikipedia.org/wiki/Feedback_arc_set )

        INPUT :

        - ``value_only`` (boolean) --
            - When set to ``True``, only the minimum
              cardinal of a minimum edge set is
              returned.

            - When set to ``False``, the ``Set`` of edges
              of a minimal edge set is returned.

        This problem is solved using Linear Programming, which certainly
        is not the best way and will have to be updated. The program solved
        is the following :

        .. MATH:
            \mbox{Minimize : }&\sum_{(u,v)\in G} b_{(u,v)}\\
            \mbox{Such that : }&\\
            &\forall v\in G, \sum_{i\in [0,\dots,n-1]}x_{v,i}=1\\
            &\forall i\in [0,\dots,n-1], \sum_{v\in G}x_{v,i}=1\\
            &\forall v\in G,\sum_{i\in [0,\dots,n-1]} ix_{v,i}=d_v\\
            &\forall (u,v)\in G, d_u-d_v+nb_{(u,v)}\geq 0\\

        An explanation :

        An acyclic digraph can be seen as a poset, and every poset has
        a linear extension. This means that in any acyclic digraph
        the vertices can be ordered with a total order `<` in such a way
        that if `(u,v)\in G`, then `u<v`.

        Thus, this linear program is built in order to assign to each vertex
        `v` an unique number `d_v\in [0,\dots,n-1]` such that if there exists
        an edge `(u,v)\in G` such that `d_v<d_u`, then the edge `(u,v)` is
        removed (`\Rightarrow x_{(u,v)}=1`).

        The number of edges removed is then minimized, which is
        the objective.

        EXAMPLE :

        If the digraph is created from a graph, and hence is symmetric
        ( if `uv` is an edge, then `vu` is an edge too ), then
        obviously the cardinality of its feedback arc set is the number
        of edges in the first graph ::

            sage: cycle=graphs.CycleGraph(5)
            sage: dcycle=DiGraph(cycle)
            sage: cycle.size()
            5
            sage: dcycle.feedback_edge_set(value_only=True)    # optional - requires GLPK or CBC
            5.0

        And in this situation, for any edge `uv` of the first graph, `uv` of `vu`
        is in the returned feedback arc set::

           sage: g = graphs.RandomGNP(5,.3)
           sage: dg = DiGraph(g)
           sage: feedback = dg.feedback_edge_set()             # optional - requires GLPK or CBC
           sage: (u,v,l) = g.edge_iterator().next()
           sage: (u,v) in feedback or (v,u) in feedback        # optional - requires GLPK or CBC
           True
        """

        from sage.numerical.mip import MixedIntegerLinearProgram

        p=MixedIntegerLinearProgram(maximization=False)

        b=p.new_variable()
        x=p.new_variable(dim=2)
        d=p.new_variable()
        n=self.order()
        N=range(n)

        # First and second constraints
        [p.add_constraint(sum([x[v][i] for i in N]),min=1,max=1) for v in self]
        [p.add_constraint(sum([x[v][i] for v in self]),min=1,max=1) for i in N]

        # Definition of d_v
        [p.add_constraint(sum([i*x[v][i] for i in N])-d[v],max=0,min=0) for v in self]

        # The removed vertices cover all the back arcs ( third condition )
        [p.add_constraint(d[u]-d[v]+n*(b[(u,v)]),min=0) for (u,v) in self.edges(labels=None)]

        p.set_binary(b)
        p.set_binary(x)

        p.set_objective(sum([b[(u,v)] for (u,v) in self.edges(labels=None)]))

        if value_only:
            return p.solve(objective_only=True)
        else:
            p.solve()

            b_sol=p.get_values(b)

            from sage.sets.set import Set
            return Set([(u,v) for (u,v) in self.edges(labels=None) if b_sol[(u,v)]==1])

    def feedback_vertex_set(self,value_only=False):
        r"""
        Computes the minimum feedback vertex set of a digraph.

        The minimum feedback vertex set of a digraph is a set of vertices
        that intersect all the circuits of the digraph.
        Equivalently, a minimum feedback vertex set of a DiGraph is a set
        `S` of vertices such that the digraph `G-S` is acyclic.

        For more informations, see
        ( http://en.wikipedia.org/wiki/Feedback_vertex_set )

        INPUT :

        - ``value_only`` (boolean) --
            - When set to ``True``, only the minimum
              cardinal of a minimum vertex set is
              returned.

            - When set to ``False``, the ``Set`` of vertices
              of a minimal feedback vertex set is returned.

        This problem is solved using Linear Programming, which certainly
        is not the best way and will have to be replaced by a better algorithm.
        The program solved is the following :

        .. MATH:
            \mbox{Minimize : }&\sum_{v\in G} b_v\\
            \mbox{Such that : }&\\
            &\forall v\in G, \sum_{i\in [0,\dots,n-1]}x_{v,i}=1\\
            &\forall i\in [0,\dots,n-1], \sum_{v\in G}x_{v,i}=1\\
            &\forall v\in G,\sum_{i\in [0,\dots,n-1]} ix_{v,i}=d_v\\
            &\forall (u,v)\in G, d_u-d_v+nb_u+nb_v\geq 0\\

        A brief explanation :

        An acyclic digraph can be seen as a poset, and every poset has
        a linear extension. This means that in any acyclic digraph
        the vertices can be ordered with a total order `<` in such a way
        that if `(u,v)\in G`, then `u<v`.
        Thus, this linear program is built in order to assign to each vertex
        `v` an unique number `d_v\in [0,\dots,n-1]` such that if there exists
        an edge `(u,v)\in G` such that `d_v<d_u`, then either `u` is removed
        (`\Rightarrow b_u=1`) or `v` is removed (`\Rightarrow b_v=1`).
        The number of vertices removed is then minimized, which is
        the objective.

        EXAMPLE:

        In a digraph built from a graph, any edge is replaced by arcs going
        in the two opposite directions, thus creating a cycle of length two.
        Hence, to remove all the cycles from the graph, each edge must see
        one of its neighbors removed : a feedback vertex set is in this
        situation a vertex cover ::

            sage: cycle=graphs.CycleGraph(5)
            sage: dcycle=DiGraph(cycle)
            sage: cycle.vertex_cover(value_only=True)         # optional - requires GLPK or CBC
            3
            sage: feedback = dcycle.feedback_vertex_set()     # optional - requires GLPK or CBC
            sage: feedback.cardinality()                      # optional - requires GLPK or CBC
            3
            sage: (u,v,l) = cycle.edge_iterator().next()
            sage: u in feedback or v in feedback              # optional - requires GLPK or CBC
            True

        For a circuit, the minimum feedback arc set is clearly `1` ::

            sage: circuit = digraphs.Circuit(5)
            sage: circuit.feedback_vertex_set(value_only=True) == 1    # optional - requires GLPK or CBC
            True
        """

        from sage.numerical.mip import MixedIntegerLinearProgram

        p=MixedIntegerLinearProgram(maximization=False)

        b=p.new_variable()
        x=p.new_variable(dim=2)
        d=p.new_variable()
        n=self.order()
        N=range(n)

        # First and second constraints
        [p.add_constraint(sum([x[v][i] for i in N]),min=1,max=1) for v in self]
        [p.add_constraint(sum([x[v][i] for v in self]),min=1,max=1) for i in N]

        # Definition of d_v
        [p.add_constraint(sum([i*x[v][i] for i in N])-d[v],max=0,min=0) for v in self]

        # The removed vertices cover all the back arcs ( third condition )
        [p.add_constraint(d[u]-d[v]+n*(b[u]+b[v]),min=0) for (u,v) in self.edges(labels=None)]

        p.set_binary(b)
        p.set_binary(x)

        p.set_objective(sum([b[v] for v in self]))

        if value_only:
            return p.solve(objective_only=True)
        else:
            p.solve()
            b_sol=p.get_values(b)

            from sage.sets.set import Set
            return Set([v for v in self if b_sol[v]==1])


    ### Construction

    def reverse(self):
        """
        Returns a copy of digraph with edges reversed in direction.

        EXAMPLES::

            sage: D = DiGraph({ 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] })
            sage: D.reverse()
            Reverse of (): Digraph on 6 vertices
        """
        H = DiGraph(multiedges=self.allows_multiple_edges(), loops=self.allows_loops())
        H.add_vertices(self)
        H.add_edges( [ (v,u,d) for (u,v,d) in self.edge_iterator() ] )
        name = self.name()
        if name is None:
            name = ''
        H.name("Reverse of (%s)"%name)
        return H

    ### Paths and cycles iterators

    def _all_paths_iterator(self, vertex, ending_vertices=None,
                            simple=False, max_length=None, trivial=False):
        r"""
        Returns an iterator over the paths of self starting with the
        given vertex.

        INPUT:

        -  ``vertex`` - the starting vertex of the paths.
        -  ``ending_vertices`` - iterable (default: None) on
           the allowed ending vertices of the paths. If None,
           then all vertices are allowed.
        -  ``simple`` - boolean (default: False). If set to True,
           then only simple paths are considered. A path is simple
           if no vertex occurs twice in it except possibly the first
           and last ones.
        -  ``max_length`` - non negative integer (default: None).
           The maximum length of the enumerated paths. If set to None,
           then all lengths are allowed.
        -  ``trivial`` - boolean (default: False). If set to True,
           then the empty paths are also enumerated.

        OUTPUT:

            iterator

        EXAMPLES::

            sage: g = DiGraph({'a' : ['a', 'b'], 'b' : ['c'], 'c' : ['d'], 'd' : ['c']}, loops=True)
            sage: pi = g._all_paths_iterator('a')
            sage: for _ in range(5): print pi.next()
            ['a', 'a']
            ['a', 'b']
            ['a', 'a', 'a']
            ['a', 'a', 'b']
            ['a', 'b', 'c']

        ::

            sage: pi = g._all_paths_iterator('b')
            sage: for _ in range(5): print pi.next()
            ['b', 'c']
            ['b', 'c', 'd']
            ['b', 'c', 'd', 'c']
            ['b', 'c', 'd', 'c', 'd']
            ['b', 'c', 'd', 'c', 'd', 'c']

        One may wish to enumerate simple paths, i.e. paths such that no vertex
        occurs twice in it except possibly the first and last one. The result
        is always finite but may be long to be computed::

            sage: pi = g._all_paths_iterator('a', simple=True)
            sage: list(pi)
            [['a', 'a'], ['a', 'b'], ['a', 'b', 'c'], ['a', 'b', 'c', 'd'], ['a', 'b', 'c', 'd', 'c']]
            sage: pi = g._all_paths_iterator('d', simple=True)
            sage: list(pi)
            [['d', 'c'], ['d', 'c', 'd']]

        It is possible to specify the allowed ending vertices::

            sage: pi = g._all_paths_iterator('a', ending_vertices=['c'])
            sage: for _ in range(5): print pi.next()
            ['a', 'b', 'c']
            ['a', 'a', 'b', 'c']
            ['a', 'a', 'a', 'b', 'c']
            ['a', 'b', 'c', 'd', 'c']
            ['a', 'a', 'a', 'a', 'b', 'c']
            sage: pi = g._all_paths_iterator('a', ending_vertices=['a', 'b'])
            sage: for _ in range(5): print pi.next()
            ['a', 'a']
            ['a', 'b']
            ['a', 'a', 'a']
            ['a', 'a', 'b']
            ['a', 'a', 'a', 'a']

        One can bound the length of the paths::

            sage: pi = g._all_paths_iterator('d', max_length=3)
            sage: list(pi)
            [['d', 'c'], ['d', 'c', 'd'], ['d', 'c', 'd', 'c']]

        Or include the trivial empty path::

            sage: pi = g._all_paths_iterator('a', max_length=3, trivial=True)
            sage: list(pi)
            [['a'], ['a', 'a'], ['a', 'b'], ['a', 'a', 'a'], ['a', 'a', 'b'], ['a', 'b', 'c'], ['a', 'a', 'a', 'a'], ['a', 'a', 'a', 'b'], ['a', 'a', 'b', 'c'], ['a', 'b', 'c', 'd']]
        """
        if ending_vertices is None:
            ending_vertices = self
        # First enumerate the empty path
        if trivial:
            yield [vertex]
        queue = [[vertex]]
        if max_length is None:
            from sage.rings.infinity import Infinity
            max_length = Infinity
        while queue:
            path = queue.pop(0)
            if len(path) > 1 and path[-1] in ending_vertices:
                yield path
            # Makes sure that the current path is not too long
            # and checks that the path is simple
            if len(path) <= max_length and (not simple or path.count(path[-1]) == 1):
                for neighbor in self.neighbor_out_iterator(path[-1]):
                    queue.append(path + [neighbor])

    def all_paths_iterator(self, starting_vertices=None, ending_vertices=None,
                           simple=False, max_length=None, trivial=False):
        r"""
        Returns an iterator over the paths of self. The paths are
        enumerated in increasing length order.

        INPUT:

        -  ``starting_vertices`` - iterable (default: None) on the
           vertices from which the paths must start. If None, then all
           vertices of the graph can be starting points.
        -  ``ending_vertices`` - iterable (default: None) on
           the allowed ending vertices of the paths. If None,
           then all vertices are allowed.
        -  ``simple`` - boolean (default: False). If set to True,
           then only simple paths are considered. A path is simple
           if no vertex occurs twice in it except if the starting
           and ending ones are equal.
        -  ``max_length`` - non negative integer (default: None).
           The maximum length of the enumerated paths. If set to None,
           then all lengths are allowed.
        -  ``trivial`` - boolean (default: False). If set to True,
           then the empty paths are also enumerated.

        OUTPUT:

            iterator

        AUTHOR:

            Alexandre Blondin Masse

        EXAMPLES::

            sage: g = DiGraph({'a' : ['a', 'b'], 'b' : ['c'], 'c' : ['d'], 'd' : ['c']}, loops=True)
            sage: pi = g.all_paths_iterator()
            sage: for _ in range(7): print pi.next()
            ['a', 'a']
            ['a', 'b']
            ['b', 'c']
            ['c', 'd']
            ['d', 'c']
            ['a', 'a', 'a']
            ['a', 'a', 'b']

        It is possible to precise the allowed starting and/or ending vertices::

            sage: pi = g.all_paths_iterator(starting_vertices=['a'])
            sage: for _ in range(5): print pi.next()
            ['a', 'a']
            ['a', 'b']
            ['a', 'a', 'a']
            ['a', 'a', 'b']
            ['a', 'b', 'c']
            sage: pi = g.all_paths_iterator(starting_vertices=['a'], ending_vertices=['b'])
            sage: for _ in range(5): print pi.next()
            ['a', 'b']
            ['a', 'a', 'b']
            ['a', 'a', 'a', 'b']
            ['a', 'a', 'a', 'a', 'b']
            ['a', 'a', 'a', 'a', 'a', 'b']

        One may prefer to enumerate only simple paths (see
        ``all_simple_paths(...)``)::

            sage: pi = g.all_paths_iterator(simple=True)
            sage: list(pi)
            [['a', 'a'], ['a', 'b'], ['b', 'c'], ['c', 'd'], ['d', 'c'], ['a', 'b', 'c'], ['b', 'c', 'd'], ['c', 'd', 'c'], ['d', 'c', 'd'], ['a', 'b', 'c', 'd'], ['b', 'c', 'd', 'c'], ['a', 'b', 'c', 'd', 'c']]

        Or simply bound the length of the enumerated paths::

            sage: pi = g.all_paths_iterator(starting_vertices=['a'], ending_vertices=['b', 'c'], max_length=6)
            sage: list(pi)
            [['a', 'b'], ['a', 'a', 'b'], ['a', 'b', 'c'], ['a', 'a', 'a', 'b'], ['a', 'a', 'b', 'c'], ['a', 'a', 'a', 'a', 'b'], ['a', 'a', 'a', 'b', 'c'], ['a', 'b', 'c', 'd', 'c'], ['a', 'a', 'a', 'a', 'a', 'b'], ['a', 'a', 'a', 'a', 'b', 'c'], ['a', 'a', 'b', 'c', 'd', 'c'], ['a', 'a', 'a', 'a', 'a', 'a', 'b'], ['a', 'a', 'a', 'a', 'a', 'b', 'c'], ['a', 'a', 'a', 'b', 'c', 'd', 'c'], ['a', 'b', 'c', 'd', 'c', 'd', 'c']]

        By default, empty paths are not enumerated, but it may be
        parametrized::

            sage: pi = g.all_paths_iterator(simple=True, trivial=True)
            sage: list(pi)
            [['a'], ['b'], ['c'], ['d'], ['a', 'a'], ['a', 'b'], ['b', 'c'], ['c', 'd'], ['d', 'c'], ['a', 'b', 'c'], ['b', 'c', 'd'], ['c', 'd', 'c'], ['d', 'c', 'd'], ['a', 'b', 'c', 'd'], ['b', 'c', 'd', 'c'], ['a', 'b', 'c', 'd', 'c']]
            sage: pi = g.all_paths_iterator(simple=True, trivial=False)
            sage: list(pi)
            [['a', 'a'], ['a', 'b'], ['b', 'c'], ['c', 'd'], ['d', 'c'], ['a', 'b', 'c'], ['b', 'c', 'd'], ['c', 'd', 'c'], ['d', 'c', 'd'], ['a', 'b', 'c', 'd'], ['b', 'c', 'd', 'c'], ['a', 'b', 'c', 'd', 'c']]
        """
        if starting_vertices is None:
            starting_vertices = self
        # We create one paths iterator per vertex
        # This is necessary if we want to iterate over paths
        # with increasing length
        vertex_iterators = dict([(v, self._all_paths_iterator(v, ending_vertices=ending_vertices, simple=simple, max_length=max_length, trivial=trivial)) for v in starting_vertices])
        paths = []
        for vi in vertex_iterators.values():
            try:
                path = vi.next()
                paths.append((len(path), path))
            except(StopIteration):
                pass
        # Since we always extract a shortest path, using a heap
        # can speed up the algorithm
        from heapq import heapify, heappop, heappush
        heapify(paths)
        while paths:
            # We choose the shortest available path
            _, shortest_path = heappop(paths)
            yield shortest_path
            # We update the path iterator to its next available path if it exists
            try:
                path = vertex_iterators[shortest_path[0]].next()
                heappush(paths, (len(path), path))
            except(StopIteration):
                pass

    def all_simple_paths(self, starting_vertices=None, ending_vertices=None,
                         max_length=None, trivial=False):
        r"""
        Returns a list of all the simple paths of self starting
        with one of the given vertices. A path is simple if no vertex
        occurs twice in it except possibly the starting and ending one.
        The paths are enumerated in increasing length order.

        INPUT:

        -  ``starting_vertices`` - list (default: None) of vertices
           from which the paths must start. If None, then all
           vertices of the graph can be starting points.
        -  ``ending_vertices`` - iterable (default: None) on
           the allowed ending vertices of the paths. If None,
           then all vertices are allowed.
        -  ``max_length`` - non negative integer (default: None).
           The maximum length of the enumerated paths. If set to None,
           then all lengths are allowed.
        -  ``trivial`` - boolean (default: False). If set to True,
           then the empty paths are also enumerated.

        OUTPUT:

            list

        .. NOTE::

            Although the number of simple paths of a finite graph
            is always finite, computing all its paths may take a very
            long time.

        EXAMPLES::

            sage: g = DiGraph({'a' : ['a', 'b'], 'b' : ['c'], 'c' : ['d'], 'd' : ['c']}, loops=True)
            sage: g.all_simple_paths()
            [['a', 'a'], ['a', 'b'], ['b', 'c'], ['c', 'd'], ['d', 'c'], ['a', 'b', 'c'], ['b', 'c', 'd'], ['c', 'd', 'c'], ['d', 'c', 'd'], ['a', 'b', 'c', 'd'], ['b', 'c', 'd', 'c'], ['a', 'b', 'c', 'd', 'c']]

        One may compute all paths having specific starting and/or
        ending vertices::

            sage: g.all_simple_paths(starting_vertices=['a'])
            [['a', 'a'], ['a', 'b'], ['a', 'b', 'c'], ['a', 'b', 'c', 'd'], ['a', 'b', 'c', 'd', 'c']]
            sage: g.all_simple_paths(starting_vertices=['a'], ending_vertices=['c'])
            [['a', 'b', 'c'], ['a', 'b', 'c', 'd', 'c']]
            sage: g.all_simple_paths(starting_vertices=['a'], ending_vertices=['b', 'c'])
            [['a', 'b'], ['a', 'b', 'c'], ['a', 'b', 'c', 'd', 'c']]

        It is also possible to bound the length of the paths::

            sage: g.all_simple_paths(max_length=2)
            [['a', 'a'], ['a', 'b'], ['b', 'c'], ['c', 'd'], ['d', 'c'], ['a', 'b', 'c'], ['b', 'c', 'd'], ['c', 'd', 'c'], ['d', 'c', 'd']]

        By default, empty paths are not enumerated, but this can
        be parametrized::

            sage: g.all_simple_paths(starting_vertices=['a'], trivial=True)
            [['a'], ['a', 'a'], ['a', 'b'], ['a', 'b', 'c'], ['a', 'b', 'c', 'd'], ['a', 'b', 'c', 'd', 'c']]
            sage: g.all_simple_paths(starting_vertices=['a'], trivial=False)
            [['a', 'a'], ['a', 'b'], ['a', 'b', 'c'], ['a', 'b', 'c', 'd'], ['a', 'b', 'c', 'd', 'c']]
        """
        return list(self.all_paths_iterator(starting_vertices=starting_vertices, ending_vertices=ending_vertices, simple=True, max_length=max_length, trivial=trivial))

    def _all_cycles_iterator_vertex(self, vertex, starting_vertices=None, simple=False,
                                    rooted=False, max_length=None, trivial=False,
                                    remove_acyclic_edges=True):
        r"""
        Returns an iterator over the cycles of self starting with the
        given vertex.

        INPUT:

        -  ``vertex`` - the starting vertex of the cycle.
        -  ``starting_vertices`` - iterable (default: None) on
           vertices from which the cycles must start. If None,
           then all vertices of the graph can be starting points.
           This argument is necessary if ``rooted`` is set to True.
        -  ``simple`` - boolean (default: False). If set to True,
           then only simple cycles are considered. A cycle is simple
           if the only vertex occuring twice in it is the starting
           and ending one.
        -  ``rooted`` - boolean (default: False). If set to False,
           then cycles differing only by their starting vertex are
           considered the same  (e.g. ``['a', 'b', 'c', 'a']`` and
           ``['b', 'c', 'a', 'b']``). Otherwise, all cycles are enumerated.
        -  ``max_length`` - non negative integer (default: None).
           The maximum length of the enumerated cycles. If set to None,
           then all lengths are allowed.
        -  ``trivial`` - boolean (default: False). If set to True,
           then the empty cycles are also enumerated.
        -  ``remove_acyclic_edges`` - boolean (default: True) which
           precises if the acyclic edges must be removed from the graph.
           Used to avoid recomputing it for each vertex.

        OUTPUT:

            iterator

        EXAMPLES::

            sage: g = DiGraph({'a' : ['a', 'b'], 'b' : ['c'], 'c' : ['d'], 'd' : ['c']}, loops=True)
            sage: it = g._all_cycles_iterator_vertex('a', simple=False, max_length=None)
            sage: for i in range(5): print it.next()
            ['a', 'a']
            ['a', 'a', 'a']
            ['a', 'a', 'a', 'a']
            ['a', 'a', 'a', 'a', 'a']
            ['a', 'a', 'a', 'a', 'a', 'a']
            sage: it = g._all_cycles_iterator_vertex('c', simple=False, max_length=None)
            sage: for i in range(5): print it.next()
            ['c', 'd', 'c']
            ['c', 'd', 'c', 'd', 'c']
            ['c', 'd', 'c', 'd', 'c', 'd', 'c']
            ['c', 'd', 'c', 'd', 'c', 'd', 'c', 'd', 'c']
            ['c', 'd', 'c', 'd', 'c', 'd', 'c', 'd', 'c', 'd', 'c']

            sage: it = g._all_cycles_iterator_vertex('d', simple=False, max_length=None)
            sage: for i in range(5): print it.next()
            ['d', 'c', 'd']
            ['d', 'c', 'd', 'c', 'd']
            ['d', 'c', 'd', 'c', 'd', 'c', 'd']
            ['d', 'c', 'd', 'c', 'd', 'c', 'd', 'c', 'd']
            ['d', 'c', 'd', 'c', 'd', 'c', 'd', 'c', 'd', 'c', 'd']

        It is possible to set a maximum length so that the number of cycles is
        finite::

            sage: it = g._all_cycles_iterator_vertex('d', simple=False, max_length=6)
            sage: list(it)
            [['d', 'c', 'd'], ['d', 'c', 'd', 'c', 'd'], ['d', 'c', 'd', 'c', 'd', 'c', 'd']]

        When ``simple`` is set to True, the number of cycles is finite since no vertex
        but the first one can occur more than once::

            sage: it = g._all_cycles_iterator_vertex('d', simple=True, max_length=None)
            sage: list(it)
            [['d', 'c', 'd']]

        By default, the empty cycle is not enumerated::

            sage: it = g._all_cycles_iterator_vertex('d', simple=True, trivial=True)
            sage: list(it)
            [['d'], ['d', 'c', 'd']]
        """
        if starting_vertices is None:
            starting_vertices = [vertex]
        # First enumerate the empty cycle
        if trivial:
            yield [vertex]
        # First we remove vertices and edges that are not part of any cycle
        if remove_acyclic_edges:
            sccs = self.strongly_connected_components()
            d = {}
            for id, component in enumerate(sccs):
                for v in component:
                    d[v] = id
            h = self.copy()
            h.delete_edges([(u,v) for (u,v) in h.edge_iterator(labels=False) if d[u] != d[v]])
        else:
            h = self
        queue = [[vertex]]
        if max_length is None:
            from sage.rings.infinity import Infinity
            max_length = Infinity
        while queue:
            path = queue.pop(0)
            # Checks if a cycle has been found
            if len(path) > 1 and path[0] == path[-1]:
                yield path
            # Makes sure that the current cycle is not too long
            # Also if a cycle has been encountered and only simple cycles are allowed,
            # Then it discards the current path
            if len(path) <= max_length and (not simple or path.count(path[-1]) == 1):
                for neighbor in h.neighbor_out_iterator(path[-1]):
                    # If cycles are not rooted, makes sure to keep only the minimum
                    # cycle according to the lexicographic order
                    if rooted or neighbor not in starting_vertices or path[0] <= neighbor:
                        queue.append(path + [neighbor])

    def all_cycles_iterator(self, starting_vertices=None, simple=False,
                            rooted=False, max_length=None, trivial=False):
        r"""
        Returns an iterator over all the cycles of self starting
        with one of the given vertices. The cycles are enumerated
        in increasing length order.

        INPUT:

        -  ``starting_vertices`` - iterable (default: None) on vertices
           from which the cycles must start. If None, then all
           vertices of the graph can be starting points.
        -  ``simple`` - boolean (default: False). If set to True,
           then only simple cycles are considered. A cycle is simple
           if the only vertex occuring twice in it is the starting
           and ending one.
        -  ``rooted`` - boolean (default: False). If set to False,
           then cycles differing only by their starting vertex are
           considered the same  (e.g. ``['a', 'b', 'c', 'a']`` and
           ``['b', 'c', 'a', 'b']``). Otherwise, all cycles are enumerated.
        -  ``max_length`` - non negative integer (default: None).
           The maximum length of the enumerated cycles. If set to None,
           then all lengths are allowed.
        -  ``trivial`` - boolean (default: False). If set to True,
           then the empty cycles are also enumerated.

        OUTPUT:

            iterator

        .. NOTE::

            See also ``all_simple_cycles(...)``.

        AUTHOR:

            Alexandre Blondin Masse

        EXAMPLES::

            sage: g = DiGraph({'a' : ['a', 'b'], 'b' : ['c'], 'c' : ['d'], 'd' : ['c']}, loops=True)
            sage: it = g.all_cycles_iterator()
            sage: for _ in range(7): print it.next()
            ['a', 'a']
            ['a', 'a', 'a']
            ['c', 'd', 'c']
            ['a', 'a', 'a', 'a']
            ['a', 'a', 'a', 'a', 'a']
            ['c', 'd', 'c', 'd', 'c']
            ['a', 'a', 'a', 'a', 'a', 'a']

        There are no cycles in the empty graph and in acyclic graphs::

            sage: g = DiGraph()
            sage: it = g.all_cycles_iterator()
            sage: list(it)
            []
            sage: g = DiGraph({0:[1]})
            sage: it = g.all_cycles_iterator()
            sage: list(it)
            []

        It is possible to restrict the starting vertices of the cycles::

            sage: g = DiGraph({'a' : ['a', 'b'], 'b' : ['c'], 'c' : ['d'], 'd' : ['c']}, loops=True)
            sage: it = g.all_cycles_iterator(starting_vertices=['b', 'c'])
            sage: for _ in range(3): print it.next()
            ['c', 'd', 'c']
            ['c', 'd', 'c', 'd', 'c']
            ['c', 'd', 'c', 'd', 'c', 'd', 'c']

        Also, one can bound the length of the cycles::

            sage: it = g.all_cycles_iterator(max_length=3)
            sage: list(it)
            [['a', 'a'], ['a', 'a', 'a'], ['c', 'd', 'c'], ['a', 'a', 'a', 'a']]

        By default, cycles differing only by their starting point are not all
        enumerated, but this may be parametrized::

            sage: it = g.all_cycles_iterator(max_length=3, rooted=False)
            sage: list(it)
            [['a', 'a'], ['a', 'a', 'a'], ['c', 'd', 'c'], ['a', 'a', 'a', 'a']]
            sage: it = g.all_cycles_iterator(max_length=3, rooted=True)
            sage: list(it)
            [['a', 'a'], ['a', 'a', 'a'], ['c', 'd', 'c'], ['d', 'c', 'd'], ['a', 'a', 'a', 'a']]

        One may prefer to enumerate simple cycles, i.e. cycles such that the only
        vertex occuring twice in it is the starting and ending one (see also
        ``all_simple_cycles(...)``::

            sage: it = g.all_cycles_iterator(simple=True)
            sage: list(it)
            [['a', 'a'], ['c', 'd', 'c']]
            sage: g = digraphs.Circuit(4)
            sage: list(g.all_cycles_iterator(simple=True))
            [[0, 1, 2, 3, 0]]
        """
        if starting_vertices is None:
            starting_vertices = self
        # Since a cycle is always included in a given strongly connected component,
        # we may remove edges from the graph
        sccs = self.strongly_connected_components()
        d = {}
        for id, component in enumerate(sccs):
            for v in component:
                d[v] = id
        h = self.copy()
        h.delete_edges([(u,v) for (u,v) in h.edge_iterator(labels=False) if d[u] != d[v]])
        # We create one cycles iterator per vertex
        # This is necessary if we want to iterate over cycles
        # with increasing length
        vertex_iterators = dict([(v, h._all_cycles_iterator_vertex(v, starting_vertices=starting_vertices, simple=simple, rooted=rooted, max_length=max_length, trivial=trivial, remove_acyclic_edges=False)) for v in starting_vertices])
        cycles = []
        for vi in vertex_iterators.values():
            try:
                cycle = vi.next()
                cycles.append((len(cycle), cycle))
            except(StopIteration):
                pass
        # Since we always extract a shortest path, using a heap
        # can speed up the algorithm
        from heapq import heapify, heappop, heappush
        heapify(cycles)
        while cycles:
            # We choose the shortest available cycle
            _, shortest_cycle = heappop(cycles)
            yield shortest_cycle
            # We update the cycle iterator to its next available cycle if it exists
            try:
                cycle = vertex_iterators[shortest_cycle[0]].next()
                heappush(cycles, (len(cycle), cycle))
            except(StopIteration):
                pass

    def all_simple_cycles(self, starting_vertices=None, rooted=False,
                          max_length=None, trivial=False):
        r"""
        Returns a list of all simple cycles of self.

        INPUT:

        -  ``starting_vertices`` - iterable (default: None) on vertices
           from which the cycles must start. If None, then all
           vertices of the graph can be starting points.
        -  ``rooted`` - boolean (default: False). If set to False,
           then equivalent cycles are merged into one single cycle
           (the one starting with minimum vertex).
           Two cycles are called equivalent if they differ only from
           their starting vertex (e.g. ``['a', 'b', 'c', 'a']`` and
           ``['b', 'c', 'a', 'b']``). Otherwise, all cycles are enumerated.
        -  ``max_length`` - non negative integer (default: None).
           The maximum length of the enumerated cycles. If set to None,
           then all lengths are allowed.
        -  ``trivial`` - boolean (default: False). If set to True,
           then the empty cycles are also enumerated.

        OUTPUT:

            list

        .. NOTE::

            Although the number of simple cycles of a finite graph
            is always finite, computing all its cycle may take a very
            long time.

        EXAMPLES::

            sage: g = DiGraph({'a' : ['a', 'b'], 'b' : ['c'], 'c' : ['d'], 'd' : ['c']}, loops=True)
            sage: g.all_simple_cycles()
            [['a', 'a'], ['c', 'd', 'c']]

        The directed version of the Petersen graph::

            sage: g = graphs.PetersenGraph().to_directed()
            sage: g.all_simple_cycles(max_length=4)
            [[0, 1, 0], [0, 4, 0], [0, 5, 0], [1, 2, 1], [1, 6, 1], [2, 3, 2], [2, 7, 2], [3, 8, 3], [3, 4, 3], [4, 9, 4], [5, 8, 5], [5, 7, 5], [6, 8, 6], [6, 9, 6], [7, 9, 7]]
            sage: g.all_simple_cycles(max_length=6)
            [[0, 1, 0], [0, 4, 0], [0, 5, 0], [1, 2, 1], [1, 6, 1], [2, 3, 2], [2, 7, 2], [3, 8, 3], [3, 4, 3], [4, 9, 4], [5, 8, 5], [5, 7, 5], [6, 8, 6], [6, 9, 6], [7, 9, 7], [0, 1, 2, 3, 4, 0], [0, 1, 2, 7, 5, 0], [0, 1, 6, 8, 5, 0], [0, 1, 6, 9, 4, 0], [0, 4, 9, 6, 1, 0], [0, 4, 9, 7, 5, 0], [0, 4, 3, 8, 5, 0], [0, 4, 3, 2, 1, 0], [0, 5, 8, 3, 4, 0], [0, 5, 8, 6, 1, 0], [0, 5, 7, 9, 4, 0], [0, 5, 7, 2, 1, 0], [1, 2, 3, 8, 6, 1], [1, 2, 7, 9, 6, 1], [1, 6, 8, 3, 2, 1], [1, 6, 9, 7, 2, 1], [2, 3, 8, 5, 7, 2], [2, 3, 4, 9, 7, 2], [2, 7, 9, 4, 3, 2], [2, 7, 5, 8, 3, 2], [3, 8, 6, 9, 4, 3], [3, 4, 9, 6, 8, 3], [5, 8, 6, 9, 7, 5], [5, 7, 9, 6, 8, 5], [0, 1, 2, 3, 8, 5, 0], [0, 1, 2, 7, 9, 4, 0], [0, 1, 6, 8, 3, 4, 0], [0, 1, 6, 9, 7, 5, 0], [0, 4, 9, 6, 8, 5, 0], [0, 4, 9, 7, 2, 1, 0], [0, 4, 3, 8, 6, 1, 0], [0, 4, 3, 2, 7, 5, 0], [0, 5, 8, 3, 2, 1, 0], [0, 5, 8, 6, 9, 4, 0], [0, 5, 7, 9, 6, 1, 0], [0, 5, 7, 2, 3, 4, 0], [1, 2, 3, 4, 9, 6, 1], [1, 2, 7, 5, 8, 6, 1], [1, 6, 8, 5, 7, 2, 1], [1, 6, 9, 4, 3, 2, 1], [2, 3, 8, 6, 9, 7, 2], [2, 7, 9, 6, 8, 3, 2], [3, 8, 5, 7, 9, 4, 3], [3, 4, 9, 7, 5, 8, 3]]

        The complete graph (without loops) on `4` vertices::

            sage: g = graphs.CompleteGraph(4).to_directed()
            sage: g.all_simple_cycles()
            [[0, 1, 0], [0, 2, 0], [0, 3, 0], [1, 2, 1], [1, 3, 1], [2, 3, 2], [0, 1, 2, 0], [0, 1, 3, 0], [0, 2, 1, 0], [0, 2, 3, 0], [0, 3, 1, 0], [0, 3, 2, 0], [1, 2, 3, 1], [1, 3, 2, 1], [0, 1, 2, 3, 0], [0, 1, 3, 2, 0], [0, 2, 1, 3, 0], [0, 2, 3, 1, 0], [0, 3, 1, 2, 0], [0, 3, 2, 1, 0]]

        If the graph contains a large number of cycles, one can bound
        the length of the cycles, or simply restrict the possible
        starting vertices of the cycles::

            sage: g = graphs.CompleteGraph(20).to_directed()
            sage: g.all_simple_cycles(max_length=2)
            [[0, 1, 0], [0, 2, 0], [0, 3, 0], [0, 4, 0], [0, 5, 0], [0, 6, 0], [0, 7, 0], [0, 8, 0], [0, 9, 0], [0, 10, 0], [0, 11, 0], [0, 12, 0], [0, 13, 0], [0, 14, 0], [0, 15, 0], [0, 16, 0], [0, 17, 0], [0, 18, 0], [0, 19, 0], [1, 2, 1], [1, 3, 1], [1, 4, 1], [1, 5, 1], [1, 6, 1], [1, 7, 1], [1, 8, 1], [1, 9, 1], [1, 10, 1], [1, 11, 1], [1, 12, 1], [1, 13, 1], [1, 14, 1], [1, 15, 1], [1, 16, 1], [1, 17, 1], [1, 18, 1], [1, 19, 1], [2, 3, 2], [2, 4, 2], [2, 5, 2], [2, 6, 2], [2, 7, 2], [2, 8, 2], [2, 9, 2], [2, 10, 2], [2, 11, 2], [2, 12, 2], [2, 13, 2], [2, 14, 2], [2, 15, 2], [2, 16, 2], [2, 17, 2], [2, 18, 2], [2, 19, 2], [3, 4, 3], [3, 5, 3], [3, 6, 3], [3, 7, 3], [3, 8, 3], [3, 9, 3], [3, 10, 3], [3, 11, 3], [3, 12, 3], [3, 13, 3], [3, 14, 3], [3, 15, 3], [3, 16, 3], [3, 17, 3], [3, 18, 3], [3, 19, 3], [4, 5, 4], [4, 6, 4], [4, 7, 4], [4, 8, 4], [4, 9, 4], [4, 10, 4], [4, 11, 4], [4, 12, 4], [4, 13, 4], [4, 14, 4], [4, 15, 4], [4, 16, 4], [4, 17, 4], [4, 18, 4], [4, 19, 4], [5, 6, 5], [5, 7, 5], [5, 8, 5], [5, 9, 5], [5, 10, 5], [5, 11, 5], [5, 12, 5], [5, 13, 5], [5, 14, 5], [5, 15, 5], [5, 16, 5], [5, 17, 5], [5, 18, 5], [5, 19, 5], [6, 7, 6], [6, 8, 6], [6, 9, 6], [6, 10, 6], [6, 11, 6], [6, 12, 6], [6, 13, 6], [6, 14, 6], [6, 15, 6], [6, 16, 6], [6, 17, 6], [6, 18, 6], [6, 19, 6], [7, 8, 7], [7, 9, 7], [7, 10, 7], [7, 11, 7], [7, 12, 7], [7, 13, 7], [7, 14, 7], [7, 15, 7], [7, 16, 7], [7, 17, 7], [7, 18, 7], [7, 19, 7], [8, 9, 8], [8, 10, 8], [8, 11, 8], [8, 12, 8], [8, 13, 8], [8, 14, 8], [8, 15, 8], [8, 16, 8], [8, 17, 8], [8, 18, 8], [8, 19, 8], [9, 10, 9], [9, 11, 9], [9, 12, 9], [9, 13, 9], [9, 14, 9], [9, 15, 9], [9, 16, 9], [9, 17, 9], [9, 18, 9], [9, 19, 9], [10, 11, 10], [10, 12, 10], [10, 13, 10], [10, 14, 10], [10, 15, 10], [10, 16, 10], [10, 17, 10], [10, 18, 10], [10, 19, 10], [11, 12, 11], [11, 13, 11], [11, 14, 11], [11, 15, 11], [11, 16, 11], [11, 17, 11], [11, 18, 11], [11, 19, 11], [12, 13, 12], [12, 14, 12], [12, 15, 12], [12, 16, 12], [12, 17, 12], [12, 18, 12], [12, 19, 12], [13, 14, 13], [13, 15, 13], [13, 16, 13], [13, 17, 13], [13, 18, 13], [13, 19, 13], [14, 15, 14], [14, 16, 14], [14, 17, 14], [14, 18, 14], [14, 19, 14], [15, 16, 15], [15, 17, 15], [15, 18, 15], [15, 19, 15], [16, 17, 16], [16, 18, 16], [16, 19, 16], [17, 18, 17], [17, 19, 17], [18, 19, 18]]
            sage: g = graphs.CompleteGraph(20).to_directed()
            sage: g.all_simple_cycles(max_length=2, starting_vertices=[0])
            [[0, 1, 0], [0, 2, 0], [0, 3, 0], [0, 4, 0], [0, 5, 0], [0, 6, 0], [0, 7, 0], [0, 8, 0], [0, 9, 0], [0, 10, 0], [0, 11, 0], [0, 12, 0], [0, 13, 0], [0, 14, 0], [0, 15, 0], [0, 16, 0], [0, 17, 0], [0, 18, 0], [0, 19, 0]]

        One may prefer to distinguish equivalent cycles having distinct
        starting vertices (compare the following examples)::

            sage: g = graphs.CompleteGraph(4).to_directed()
            sage: g.all_simple_cycles(max_length=2, rooted=False)
            [[0, 1, 0], [0, 2, 0], [0, 3, 0], [1, 2, 1], [1, 3, 1], [2, 3, 2]]
            sage: g.all_simple_cycles(max_length=2, rooted=True)
            [[0, 1, 0], [0, 2, 0], [0, 3, 0], [1, 0, 1], [1, 2, 1], [1, 3, 1], [2, 0, 2], [2, 1, 2], [2, 3, 2], [3, 0, 3], [3, 1, 3], [3, 2, 3]]
        """
        return list(self.all_cycles_iterator(starting_vertices=starting_vertices, simple=True, rooted=rooted, max_length=max_length, trivial=trivial))

    ### Directed Acyclic Graphs (DAGs)

    def topological_sort(self):
        """
        Returns a topological sort of the digraph if it is acyclic, and
        raises a TypeError if the digraph contains a directed cycle.

        A topological sort is an ordering of the vertices of the digraph
        such that each vertex comes before all of its successors. That is,
        if u comes before v in the sort, then there may be a directed path
        from u to v, but there will be no directed path from v to u.

        EXAMPLES::

            sage: D = DiGraph({ 0:[1,2,3], 4:[2,5], 1:[8], 2:[7], 3:[7], 5:[6,7], 7:[8], 6:[9], 8:[10], 9:[10] })
            sage: D.plot(layout='circular').show()
            sage: D.topological_sort()
            [4, 5, 6, 9, 0, 1, 2, 3, 7, 8, 10]

        ::

            sage: D.add_edge(9,7)
            sage: D.topological_sort()
            [4, 5, 6, 9, 0, 1, 2, 3, 7, 8, 10]

        ::

            sage: D.add_edge(7,4)
            sage: D.topological_sort()
            Traceback (most recent call last):
            ...
            TypeError: Digraph is not acyclic-- there is no topological sort.

        .. note::

           There is a recursive version of this in NetworkX, but it has
           problems::

              sage: import networkx
              sage: D = DiGraph({ 0:[1,2,3], 4:[2,5], 1:[8], 2:[7], 3:[7], 5:[6,7], 7:[8], 6:[9], 8:[10], 9:[10] })
              sage: N = D.networkx_graph()
              sage: networkx.topological_sort(N)
              [4, 5, 6, 9, 0, 1, 2, 3, 7, 8, 10]
              sage: networkx.topological_sort_recursive(N) is None
              True
        """
        import networkx
        S = networkx.topological_sort(self.networkx_graph(copy=False))
        if S is None:
            raise TypeError('Digraph is not acyclic-- there is no topological sort.')
        else:
            return S

    def topological_sort_generator(self):
        """
        Returns a list of all topological sorts of the digraph if it is
        acyclic, and raises a TypeError if the digraph contains a directed
        cycle.

        A topological sort is an ordering of the vertices of the digraph
        such that each vertex comes before all of its successors. That is,
        if u comes before v in the sort, then there may be a directed path
        from u to v, but there will be no directed path from v to u. See
        also Graph.topological_sort().

        AUTHORS:

        - Mike Hansen - original implementation

        - Robert L. Miller: wrapping, documentation

        REFERENCE:

        - [1] Pruesse, Gara and Ruskey, Frank. Generating Linear
          Extensions Fast. SIAM J. Comput., Vol. 23 (1994), no. 2, pp.
          373-386.

        EXAMPLES::

            sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: D.plot(layout='circular').show()
            sage: D.topological_sort_generator()
            [[0, 1, 2, 3, 4], [0, 1, 2, 4, 3], [0, 2, 1, 3, 4], [0, 2, 1, 4, 3], [0, 2, 4, 1, 3]]

        ::

            sage: for sort in D.topological_sort_generator():
            ...       for edge in D.edge_iterator():
            ...           u,v,l = edge
            ...           if sort.index(u) > sort.index(v):
            ...               print "This should never happen."
        """
        from sage.graphs.linearextensions import LinearExtensions
        try:
            return LinearExtensions(self).list()
        except TypeError:
            raise TypeError('Digraph is not acyclic-- there is no topological sort (or there was an error in sage/graphs/linearextensions.py).')

    ### Visualization

    def graphviz_string(self):
       r"""
       Returns a representation in the DOT language, ready to render in
       graphviz.

       REFERENCES:

       - http://www.graphviz.org/doc/info/lang.html

       EXAMPLES::

           sage: G = DiGraph({0:{1:None,2:None}, 1:{2:None}, 2:{3:'foo'}, 3:{}} ,sparse=True)
           sage: s = G.graphviz_string(); s
           'digraph {\n"0";"1";"2";"3";\n"0"->"1";"0"->"2";"1"->"2";"2"->"3"[label="foo"];\n}'
       """
       return self._graphviz_string_helper("digraph", "->") # edge_string is "->" for directed graphs


    def strongly_connected_components(self):
        """
        Returns a list of lists of vertices, each list representing a
        strongly connected component.

        EXAMPLES::

            sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: D.connected_components()
            [[0, 1, 2, 3], [4, 5, 6]]
            sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: D.strongly_connected_components()
            [[0], [1], [2], [3], [4], [5], [6]]
            sage: D.add_edge([2,0])
            sage: D.strongly_connected_components()
            [[0, 1, 2], [3], [4], [5], [6]]
        """

        try:
            vertices = set(self.vertices())
            scc = []
            while vertices:
                tmp = self.strongly_connected_component_containing_vertex(vertices.__iter__().next())
                vertices.difference_update(set(tmp))
                scc.append(tmp)
            return scc

        except ValueError:
            import networkx
            return networkx.strongly_connected_components(self.networkx_graph(copy=False))


    def strongly_connected_component_containing_vertex(self, v):
        """
        Returns the set of vertices whose strongly connected component is the same as v's.

        INPUT:

        - ``v`` -- a vertex

        EXAMPLE:

        In the symmetric digraph of a graph, the strongly connected components are the connected
        components::

            sage: g = graphs.PetersenGraph()
            sage: d = DiGraph(g)
            sage: d.strongly_connected_component_containing_vertex(0)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        """

        if self.order()==1:
            return [v]

        try:
            return self._backend.strongly_connected_component_containing_vertex(v)

        except AttributeError:
            raise ValueError("This function is only defined for C graphs.")

    def strongly_connected_components_subgraphs(self):
        r"""
        Returns the strongly connected components as a list of subgraphs.

        EXAMPLE:

        In the symmetric digraph of a graph, the strongly connected components are the connected
        components::

            sage: g = graphs.PetersenGraph()
            sage: d = DiGraph(g)
            sage: d.strongly_connected_components_subgraphs()
            [Subgraph of (Petersen graph): Digraph on 10 vertices]

        """
        return map(self.subgraph, self.strongly_connected_components())


    def strongly_connected_components_digraph(self):
        r"""
        Returns the digraph of the strongly connected components

        The digraph of the strongly connected components of a graph `G` has
        a vertex per strongly connected component included in `G`. There
        is an edge from a component `C_1` to a component `C_2` if there is
        an edge from one to the other in `G`.

        EXAMPLE:

        Such a digraph is always acyclic ::

            sage: g = digraphs.RandomDirectedGNP(15,.1)
            sage: scc_digraph = g.strongly_connected_components_digraph()
            sage: scc_digraph.is_directed_acyclic()
            True

        The vertices of the digraph of strongly connected components are
        exactly the strongly connected components::

            sage: g = digraphs.ButterflyGraph(2)
            sage: scc_digraph = g.strongly_connected_components_digraph()
            sage: g.is_directed_acyclic()
            True
            sage: all([ Set(scc) in scc_digraph.vertices() for scc in g.strongly_connected_components()])
            True
        """

        from sage.sets.set import Set

        scc = self.strongly_connected_components()
        scc_set = map(Set,scc)

        d = {}
        for i,c in enumerate(scc):
            for v in c:
                d[v] = i

        g = DiGraph()
        g.add_vertices(scc_set)
        g.allow_multiple_edges(False)

        g.add_edges( (scc_set[d[u]], scc_set[d[u]]) for u,v in self.edges(labels=False) )

        return g



    def is_strongly_connected(self):
        r"""
        Returns whether the current ``DiGraph`` is strongly connected.

        EXAMPLE:

        The circuit is obviously strongly connected ::

            sage: g = digraphs.Circuit(5)
            sage: g.is_strongly_connected()
            True

        But a transitive triangle is not::

            sage: g = DiGraph({ 0 : [1,2], 1 : [2]})
            sage: g.is_strongly_connected()
            False
        """
        if self.order()==1:
            return True

        try:
            return self._backend.is_strongly_connected()

        except AttributeError:
            return len(self.strongly_connected_components()) == 1



