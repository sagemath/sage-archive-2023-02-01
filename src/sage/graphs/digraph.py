r"""
Directed graphs

This module implements functions and operations involving directed
graphs. Here is what they can do

**Graph basic operations:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~DiGraph.layout_acyclic_dummy` | Computes a (dummy) ranked layout so that all edges point upward.
    :meth:`~DiGraph.layout_acyclic` | Computes a ranked layout so that all edges point upward.
    :meth:`~DiGraph.reverse` | Returns a copy of digraph with edges reversed in direction.
    :meth:`~DiGraph.reverse_edge` | Reverses an edge.
    :meth:`~DiGraph.reverse_edges` | Reverses a list of edges.
    :meth:`~DiGraph.out_degree_sequence` | Return the outdegree sequence.
    :meth:`~DiGraph.out_degree_iterator` | Same as degree_iterator, but for out degree.
    :meth:`~DiGraph.out_degree` | Same as degree, but for out degree.
    :meth:`~DiGraph.in_degree_sequence` | Return the indegree sequence of this digraph.
    :meth:`~DiGraph.in_degree_iterator` | Same as degree_iterator, but for in degree.
    :meth:`~DiGraph.in_degree` | Same as degree, but for in-degree.
    :meth:`~DiGraph.neighbors_out` | Returns the list of the out-neighbors of a given vertex.
    :meth:`~DiGraph.neighbor_out_iterator` | Returns an iterator over the out-neighbors of a given vertex.
    :meth:`~DiGraph.neighbors_in` | Returns the list of the in-neighbors of a given vertex.
    :meth:`~DiGraph.neighbor_in_iterator` | Returns an iterator over the in-neighbors of vertex.
    :meth:`~DiGraph.outgoing_edges` | Returns a list of edges departing from vertices.
    :meth:`~DiGraph.outgoing_edge_iterator` | Return an iterator over all departing edges from vertices
    :meth:`~DiGraph.incoming_edges` | Returns a list of edges arriving at vertices.
    :meth:`~DiGraph.incoming_edge_iterator` | Return an iterator over all arriving edges from vertices
    :meth:`~DiGraph.sources` | Returns the list of all sources (vertices without incoming edges) of this digraph.
    :meth:`~DiGraph.sinks` | Returns the list of all sinks (vertices without outoing edges) of this digraph.
    :meth:`~DiGraph.to_undirected` | Returns an undirected version of the graph.
    :meth:`~DiGraph.to_directed` | Since the graph is already directed, simply returns a copy of itself.
    :meth:`~DiGraph.is_directed` | Since digraph is directed, returns True.
    :meth:`~DiGraph.dig6_string` | Returns the dig6 representation of the digraph as an ASCII string.

**Paths and cycles:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~DiGraph.all_paths_iterator` | Returns an iterator over the paths of self. The paths are
    :meth:`~DiGraph.all_simple_paths` | Returns a list of all the simple paths of self starting
    :meth:`~DiGraph.all_cycles_iterator` | Returns an iterator over all the cycles of self starting
    :meth:`~DiGraph.all_simple_cycles` | Returns a list of all simple cycles of self.

**Representation theory:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~Digraph.path_semigroup` | Returns the (partial) semigroup formed by the paths of the digraph.

**Connectivity:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~DiGraph.is_strongly_connected` | Returns whether the current ``DiGraph`` is strongly connected.
    :meth:`~DiGraph.strongly_connected_components_digraph` | Returns the digraph of the strongly connected components
    :meth:`~DiGraph.strongly_connected_components_subgraphs` | Returns the strongly connected components as a list of subgraphs.
    :meth:`~DiGraph.strongly_connected_component_containing_vertex` | Returns the strongly connected component containing a given vertex
    :meth:`~DiGraph.strongly_connected_components` | Returns the list of strongly connected components.


**Acyclicity:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~DiGraph.is_directed_acyclic` | Returns whether the digraph is acyclic or not.
    :meth:`~DiGraph.is_transitive` | Returns whether the digraph is transitive or not.
    :meth:`~DiGraph.is_aperiodic` | Returns whether the digraph is aperiodic or not.
    :meth:`~DiGraph.level_sets` | Returns the level set decomposition of the digraph.
    :meth:`~DiGraph.topological_sort_generator` | Returns a list of all topological sorts of the digraph if it is acyclic
    :meth:`~DiGraph.topological_sort` | Returns a topological sort of the digraph if it is acyclic

**Hard stuff:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~DiGraph.feedback_edge_set` | Computes the minimum feedback edge (arc) set of a digraph

Methods
-------
"""

from sage.rings.integer import Integer
from sage.misc.superseded import deprecated_function_alias
from sage.misc.superseded import deprecation
import sage.graphs.generic_graph_pyx as generic_graph_pyx
from sage.graphs.generic_graph import GenericGraph
from sage.graphs.dot2tex_utils import have_dot2tex

class DiGraph(GenericGraph):
    """
    Directed graph.

    A digraph or directed graph is a set of vertices connected by oriented
    edges. For more information, see the
    `Wikipedia article on digraphs
    <http://en.wikipedia.org/wiki/Digraph_%28mathematics%29>`_.

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

    -  ``data`` -  can be any of the following (see the ``format`` keyword):

       #.  A dictionary of dictionaries

       #.  A dictionary of lists

       #.  A numpy matrix or ndarray

       #.  A Sage adjacency matrix or incidence matrix

       #.  A pygraphviz graph

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

       -  ``NX`` - data must be a NetworkX DiGraph.

           .. NOTE::

               As Sage's default edge labels is ``None`` while NetworkX uses
               ``{}``, the ``{}`` labels of a NetworkX digraph are automatically
               set to ``None`` when it is converted to a Sage graph. This
               behaviour can be overruled by setting the keyword
               ``convert_empty_dict_labels_to_None`` to ``False`` (it is
               ``True`` by default).

    -  ``boundary`` - a list of boundary vertices, if
       empty, digraph is considered as a 'graph without boundary'

    -  ``implementation`` - what to use as a backend for
       the graph. Currently, the options are either 'networkx' or
       'c_graph'

    - ``sparse`` (boolean) -- ``sparse=True`` is an alias for
      ``data_structure="sparse"``, and ``sparse=False`` is an alias for
      ``data_structure="dense"``.

    -  ``data_structure`` -- one of the following

       * ``"dense"`` -- selects the :mod:`~sage.graphs.base.dense_graph`
         backend.

       * ``"sparse"`` -- selects the :mod:`~sage.graphs.base.sparse_graph`
         backend.

       * ``"static_sparse"`` -- selects the
         :mod:`~sage.graphs.base.static_sparse_backend` (this backend is faster
         than the sparse backend and smaller in memory, and it is immutable, so
         that the resulting graphs can be used as dictionary keys).

       *Only available when* ``implementation == 'c_graph'``

    - ``immutable`` (boolean) -- whether to create a immutable digraph. Note
      that ``immutable=True`` is actually a shortcut for
      ``data_structure='static_sparse'``. Set to ``False`` by default, only
      available when ``implementation='c_graph'``

    -  ``vertex_labels`` - only for implementation == 'c_graph'.
       Whether to allow any object as a vertex (slower), or
       only the integers 0, ..., n-1, where n is the number of vertices.

    -  ``convert_empty_dict_labels_to_None`` - this arguments sets
       the default edge labels used by NetworkX (empty dictionaries)
       to be replaced by None, the default Sage edge label. It is
       set to ``True`` iff a NetworkX graph is on the input.

    EXAMPLES:

    #. A dictionary of dictionaries::

            sage: g = DiGraph({0:{1:'x',2:'z',3:'a'}, 2:{5:'out'}}); g
            Digraph on 5 vertices

       The labels ('x', 'z', 'a', 'out') are labels for edges. For
       example, 'out' is the label for the edge from 2 to 5. Labels can be
       used as weights, if all the labels share some common parent.

    #. A dictionary of lists (or iterables)::

            sage: g = DiGraph({0:[1,2,3], 2:[4]}); g
            Digraph on 5 vertices
            sage: g = DiGraph({0:(1,2,3), 2:(4,)}); g
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

    TESTS::

        sage: DiGraph({0:[1,2,3], 2:[4]}).edges()
        [(0, 1, None), (0, 2, None), (0, 3, None), (2, 4, None)]
        sage: DiGraph({0:(1,2,3), 2:(4,)}).edges()
        [(0, 1, None), (0, 2, None), (0, 3, None), (2, 4, None)]
        sage: DiGraph({0:Set([1,2,3]), 2:Set([4])}).edges()
        [(0, 1, None), (0, 2, None), (0, 3, None), (2, 4, None)]

    Get rid of mutable default argument for `boundary` (:trac:`14794`)::

        sage: D = DiGraph(boundary=None)
        sage: D._boundary
        []

    Demonstrate that digraphs using the static backend are equal to mutable
    graphs but can be used as dictionary keys::

        sage: import networkx
        sage: g = networkx.DiGraph({0:[1,2,3], 2:[4]})
        sage: G = DiGraph(g, implementation='networkx')
        sage: G_imm = DiGraph(G, data_structure="static_sparse")
        sage: H_imm = DiGraph(G, data_structure="static_sparse")
        sage: H_imm is G_imm
        False
        sage: H_imm == G_imm == G
        True
        sage: {G_imm:1}[H_imm]
        1
        sage: {G_imm:1}[G]
        Traceback (most recent call last):
        ...
        TypeError: This graph is mutable, and thus not hashable. Create an
        immutable copy by `g.copy(immutable=True)`

    The error message states that one can also create immutable graphs by
    specifying the ``immutable`` optional argument (not only by
    ``data_structure='static_sparse'`` as above)::

        sage: J_imm = DiGraph(G, immutable=True)
        sage: J_imm == G_imm
        True
        sage: type(J_imm._backend) == type(G_imm._backend)
        True

    """
    _directed = True

    def __init__(self, data=None, pos=None, loops=None, format=None,
                 boundary=None, weighted=None, implementation='c_graph',
                 data_structure="sparse", vertex_labels=True, name=None,
                 multiedges=None, convert_empty_dict_labels_to_None=None,
                 sparse=True, immutable=False):
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

        Invalid sequence of edges given as an input (they do not all
        have the same length)::

            sage: g = DiGraph([(1,2),(2,3),(2,3,4)])
            Traceback (most recent call last):
            ...
            ValueError: Edges input must all follow the same format.


        Two different labels given for the same edge in a graph
        without multiple edges::

            sage: g = Graph([(1,2,3),(1,2,4)], multiedges = False)
            Traceback (most recent call last):
            ...
            ValueError: Two different labels given for the same edge in a graph without multiple edges.

        Detection of multiple edges::

            sage: DiGraph([(1, 2, 0), (1,2,1)])
            Multi-digraph on 2 vertices
            sage: DiGraph([(1, 2, 0)])
            Digraph on 2 vertices

        An empty list or dictionary defines a simple graph (trac #10441 and #12910)::

            sage: DiGraph([])
            Digraph on 0 vertices
            sage: DiGraph({})
            Digraph on 0 vertices
            sage: # not "Multi-digraph on 0 vertices"

        Problem with weighted adjacency matrix (:trac:`13919`)::

            sage: B = {0:{1:2,2:5,3:4},1:{2:2,4:7},2:{3:1,4:4,5:3},3:{5:4},4:{5:1,6:5},5:{4:1,6:7,5:1}}
            sage: grafo3 = DiGraph(B,weighted=True)
            sage: matad = grafo3.weighted_adjacency_matrix()
            sage: grafo4 = DiGraph(matad,format = "adjacency_matrix", weighted=True)
            sage: grafo4.shortest_path(0,6,by_weight=True)
            [0, 1, 2, 5, 4, 6]

        Building a DiGraph with ``immutable=False`` returns a mutable graph::

            sage: g = graphs.PetersenGraph()
            sage: g = DiGraph(g.edges(),immutable=False)
            sage: g.add_edge("Hey", "Heyyyyyyy")
            sage: {g:1}[g]
            Traceback (most recent call last):
            ...
            TypeError: This graph is mutable, and thus not hashable. Create an immutable copy by `g.copy(immutable=True)`
            sage: copy(g) is g
            False
            sage: {g.copy(immutable=True):1}[g.copy(immutable=True)]
            1

        But building it with ``immutable=True`` returns an immutable graph::

            sage: g = DiGraph(graphs.PetersenGraph(), immutable=True)
            sage: g.add_edge("Hey", "Heyyyyyyy")
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: {g:1}[g]
            1
            sage: copy(g) is g
            True

        Check the error when multiple edges are sent but ``multiple_edges`` is
        set to ``False`` (:trac:`16215`)::

            sage: DiGraph([(0,1),(1,0),(0,1)], multiedges=False)
            Traceback (most recent call last):
            ...
            ValueError: Non-multidigraph got several edges (0,1)
        """
        msg = ''
        GenericGraph.__init__(self)
        from sage.structure.element import is_Matrix
        from sage.misc.misc import uniq

        if sparse == False:
            if data_structure != "sparse":
                raise ValueError("The 'sparse' argument is an alias for "
                                 "'data_structure'. Please do not define both.")
            data_structure = "dense"

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
                if isinstance(data[keys[0]], dict):
                    format = 'dict_of_dicts'
                else:
                    format = 'dict_of_lists'
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

        # Input is a list of edges
        if format is None and isinstance(data, list):

            # If we are given a list (we assume it is a list of
            # edges), we convert the problem to a dict_of_dicts or a
            # dict_of_lists

            edges = data

            # Along the process, we are assuming all edges have the
            # same length. If it is not true, a ValueError will occur
            try:

                # The list is empty
                if not data:
                    data = {}
                    format = 'dict_of_dicts'

                # The edges are not labelled
                elif len(data[0]) == 2:
                    data = {}
                    for u,v in edges:
                        if not u in data:
                            data[u] = []

                        data[u].append(v)

                    format = 'dict_of_lists'

                # The edges are labelled
                elif len(data[0]) == 3:
                    data = {}
                    for u,v,l in edges:

                        if not u in data:
                            data[u] = {}

                        # Now the key exists, and is a dictionary ...

                        # If we notice for the first time that there
                        # are multiple edges, we update the whole
                        # dictionary so that data[u][v] is a list

                        if (multiedges is None and
                            (v in data[u])):
                            multiedges = True
                            for uu, dd in data.iteritems():
                                for vv, ddd in dd.iteritems():
                                    dd[vv] = [ddd]

                        # If multiedges is set to False while the same
                        # edge is present in the list with different
                        # values of its label
                        elif (multiedges is False and
                            (v in data[u] and data[u][v] != l)):
                                raise ValueError("MULTIEDGE")

                        # Now we are behaving as if multiedges == None
                        # means multiedges = False. If something bad
                        # happens later, the whole dictionary will be
                        # updated anyway

                        if multiedges is True:
                            if v not in data[u]:
                                data[u][v] = []

                            data[u][v].append(l)

                        else:
                            data[u][v] = l

                    format = 'dict_of_dicts'

            except ValueError as ve:
                if str(ve) == "MULTIEDGE":
                    raise ValueError("Two different labels given for the same edge in a graph without multiple edges.")
                else:
                    raise ValueError("Edges input must all follow the same format.")

        if format is None:
            import networkx
            data = networkx.MultiDiGraph(data)
            format = 'NX'

        # At this point, format has been set.

        # adjust for empty dicts instead of None in NetworkX default edge labels
        if convert_empty_dict_labels_to_None is None:
            convert_empty_dict_labels_to_None = (format == 'NX')

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
                except Exception:
                    if weighted is False:
                        raise ValueError("Non-weighted digraph's"+
                        " adjacency matrix must have only nonnegative"+
                        " integer entries")
                    weighted = True
                    if multiedges is None: multiedges = False
                    break

            if weighted is None:
                weighted = False

            if multiedges is None:
                multiedges = ((not weighted) and sorted(entries) != [0,1])

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
                    L = sorted(uniq(c.list()))
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
            if multiedges is None and len(data) > 0:
                multiedges = True
            num_verts = len(verts)
        elif format == 'dict_of_lists':
            # convert to a dict of lists if not already one
            if not all(isinstance(data[u], list) for u in data):
                data = dict([u, list(data[u])] for u in data)
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
                        v = (v for v in data[u] if data[u].count(v) > 1).next()
                        raise ValueError("Non-multidigraph got several edges (%s,%s)"%(u,v))
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
                if data_structure == "dense":
                    raise RuntimeError("Multiedge and weighted c_graphs must be sparse.")

            if immutable:
                data_structure = 'static_sparse'

            # If the data structure is static_sparse, we first build a graph
            # using the sparse data structure, then reencode the resulting graph
            # as a static sparse graph.
            from sage.graphs.base.sparse_graph import SparseGraphBackend
            from sage.graphs.base.dense_graph import DenseGraphBackend
            if data_structure in ["sparse", "static_sparse"]:
                CGB = SparseGraphBackend
            elif data_structure == "dense":
                 CGB = DenseGraphBackend
            else:
                raise ValueError("data_structure must be equal to 'sparse', "
                                 "'static_sparse' or 'dense'")
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
        else:
            raise NotImplementedError("Supported implementations: networkx, c_graph.")

        if format == 'dig6':
            k = 0
            for i in xrange(n):
                for j in xrange(n):
                    if m[k] == '1':
                        self._backend.add_edge(i, j, None, True)
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
                        self._backend.add_edge(uu,vv,None,True)
        elif format == 'dict_of_dicts':
            if convert_empty_dict_labels_to_None:
                for u in data:
                    for v in data[u]:
                        if multiedges:
                            for l in data[u][v]:
                                self._backend.add_edge(u,v,l,True)
                        else:
                            self._backend.add_edge(u,v,data[u][v] if data[u][v] != {} else None,True)
            else:
                for u in data:
                    for v in data[u]:
                        if multiedges:
                            for l in data[u][v]:
                                self._backend.add_edge(u,v,l,True)
                        else:
                            self._backend.add_edge(u,v,data[u][v],True)
        elif format == 'dict_of_lists':
            for u in data:
                for v in data[u]:
                    self._backend.add_edge(u,v,None,True)
        else:
            assert format == 'int'
        self._pos = pos
        self._boundary = boundary if boundary is not None else []
        if format != 'DiGraph' or name is not None:
            self.name(name)

        if data_structure == "static_sparse":
            from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            ib = StaticSparseBackend(self, loops = loops, multiedges = multiedges)
            self._backend = ib
            self._immutable = True

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

    def is_directed_acyclic(self, certificate = False):
        """
        Returns whether the digraph is acyclic or not.

        A directed graph is acyclic if for any vertex `v`, there is no directed
        path that starts and ends at `v`. Every directed acyclic graph (DAG)
        corresponds to a partial ordering of its vertices, however multiple dags
        may lead to the same partial ordering.

        INPUT:

        - ``certificate`` -- whether to return a certificate (``False`` by
          default).

        OUTPUT:

            * When ``certificate=False``, returns a boolean value.

            * When ``certificate=True`` :

                * If the graph is acyclic, returns a pair ``(True, ordering)``
                  where ``ordering`` is a list of the vertices such that ``u``
                  appears before ``v`` in ``ordering`` if ``u, v`` is an edge.

                * Else, returns a pair ``(False, cycle)`` where ``cycle`` is a
                  list of vertices representing a circuit in the graph.

        EXAMPLES:

        At first, the following graph is acyclic::

            sage: D = DiGraph({ 0:[1,2,3], 4:[2,5], 1:[8], 2:[7], 3:[7], 5:[6,7], 7:[8], 6:[9], 8:[10], 9:[10] })
            sage: D.plot(layout='circular').show()
            sage: D.is_directed_acyclic()
            True

        Adding an edge from `9` to `7` does not change it::

            sage: D.add_edge(9,7)
            sage: D.is_directed_acyclic()
            True

        We can obtain as a proof an ordering of the vertices such that `u`
        appears before `v` if `uv` is an edge of the graph::

            sage: D.is_directed_acyclic(certificate = True)
            (True, [4, 5, 6, 9, 0, 1, 2, 3, 7, 8, 10])

        Adding an edge from 7 to 4, though, makes a difference::

            sage: D.add_edge(7,4)
            sage: D.is_directed_acyclic()
            False

        Indeed, it creates a circuit `7, 4, 5`::

            sage: D.is_directed_acyclic(certificate = True)
            (False, [7, 4, 5])

        Checking acyclic graphs are indeed acyclic ::

            sage: def random_acyclic(n, p):
            ...    g = graphs.RandomGNP(n, p)
            ...    h = DiGraph()
            ...    h.add_edges([ ((u,v) if u<v else (v,u)) for u,v,_ in g.edges() ])
            ...    return h
            ...
            sage: all( random_acyclic(100, .2).is_directed_acyclic()    # long time
            ...        for i in range(50))                              # long time
            True

        TESTS:

        What about loops?::

            sage: g = digraphs.ButterflyGraph(3)
            sage: g.allow_loops(True)
            sage: g.is_directed_acyclic()
            True
            sage: g.add_edge(0,0)
            sage: g.is_directed_acyclic()
            False
        """
        return self._backend.is_directed_acyclic(certificate = certificate)

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

    def to_undirected(self, implementation='c_graph', data_structure=None,
                      sparse=None):
        """
        Returns an undirected version of the graph. Every directed edge
        becomes an edge.

        INPUT:

         - ``implementation`` - string (default: 'networkx') the
           implementation goes here.  Current options are only
           'networkx' or 'c_graph'.

         - ``data_structure`` -- one of ``"sparse"``, ``"static_sparse"``, or
           ``"dense"``. See the documentation of :class:`Graph` or
           :class:`DiGraph`.

         - ``sparse`` (boolean) -- ``sparse=True`` is an alias for
           ``data_structure="sparse"``, and ``sparse=False`` is an alias for
           ``data_structure="dense"``.

        EXAMPLES::

            sage: D = DiGraph({0:[1,2],1:[0]})
            sage: G = D.to_undirected()
            sage: D.edges(labels=False)
            [(0, 1), (0, 2), (1, 0)]
            sage: G.edges(labels=False)
            [(0, 1), (0, 2)]
        """
        if sparse is not None:
            deprecation(14806,"The 'sparse' keyword has been deprecated, and "
                        "is now replaced by 'data_structure' which has a different "
                        "meaning. Please consult the documentation.")
            data_structure = "sparse" if sparse else "dense"

        if data_structure is None:
            from sage.graphs.base.dense_graph import DenseGraphBackend
            from sage.graphs.base.sparse_graph import SparseGraphBackend
            if isinstance(self._backend, DenseGraphBackend):
                data_structure = "dense"
            elif isinstance(self._backend, SparseGraphBackend):
                data_structure = "sparse"
            else:
                data_structure = "static_sparse"
        from sage.graphs.all import Graph
        G = Graph(name=self.name(), pos=self._pos, boundary=self._boundary,
                  multiedges=self.allows_multiple_edges(), loops=self.allows_loops(),
                  implementation=implementation, data_structure=data_structure)
        G.name(self.name())
        G.add_vertices(self.vertex_iterator())
        G.add_edges(self.edge_iterator())
        if hasattr(self, '_embedding'):
            from copy import copy
            G._embedding = copy(self._embedding)
        G._weighted = self._weighted
        return G

    ### Edge Handlers

    def incoming_edge_iterator(self, vertices, labels=True):
        """
        Return an iterator over all arriving edges from vertices.

        INPUT:

        - ``vertices`` -- a vertex or a list of vertices

        - ``labels`` (boolean) -- whether to return edges as pairs of vertices,
          or as triples containing the labels.

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

    def incoming_edges(self, vertices, labels=True):
        """
        Returns a list of edges arriving at vertices.

        INPUT:

        - ``vertices`` -- a vertex or a list of vertices

        - ``labels`` (boolean) -- whether to return edges as pairs of vertices,
          or as triples containing the labels.

        EXAMPLES::

            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.incoming_edges([0])
            [(1, 0, None), (4, 0, None)]
        """
        return list(self.incoming_edge_iterator(vertices, labels=labels))

    def outgoing_edge_iterator(self, vertices, labels=True):
        """
        Return an iterator over all departing edges from vertices.

        INPUT:

        - ``vertices`` -- a vertex or a list of vertices

        - ``labels`` (boolean) -- whether to return edges as pairs of vertices,
          or as triples containing the labels.

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

    def outgoing_edges(self, vertices, labels=True):
        """
        Returns a list of edges departing from vertices.

        INPUT:

        - ``vertices`` -- a vertex or a list of vertices

        - ``labels`` (boolean) -- whether to return edges as pairs of vertices,
          or as triples containing the labels.

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

    predecessor_iterator = deprecated_function_alias(7634, neighbor_in_iterator)

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

    predecessors = deprecated_function_alias(7634, neighbors_in)

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

    successor_iterator = deprecated_function_alias(7634, neighbor_out_iterator)

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

    successors = deprecated_function_alias(7634, neighbors_out)

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
        Return the indegree sequence.

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
            sage: D.out_degree(2)
            1
        """
        if vertices in self:

            try:
                return self._backend.out_degree(vertices)

            except AttributeError:
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

    def sources(self):
        r"""
        Returns a list of sources of the digraph.

        OUTPUT:

        - list, the vertices of the digraph that have no edges going into them

        EXAMPLES::

            sage: G = DiGraph({1:{3:['a']}, 2:{3:['b']}})
            sage: G.sources()
            [1, 2]
            sage: T = DiGraph({1:{}})
            sage: T.sources()
            [1]
        """
        return [x for x in self if self.in_degree(x)==0]

    def sinks(self):
        """
        Returns a list of sinks of the digraph.

        OUTPUT:

        - list, the vertices of the digraph that have no edges beginning at them

        EXAMPLES::

            sage: G = DiGraph({1:{3:['a']}, 2:{3:['b']}})
            sage: G.sinks()
            [3]
            sage: T = DiGraph({1:{}})
            sage: T.sinks()
            [1]
        """
        return [x for x in self if self.out_degree(x)==0]


    def feedback_edge_set(self, constraint_generation= True, value_only=False, solver=None, verbose=0):
        r"""
        Computes the minimum feedback edge set of a digraph (also called
        feedback arc set).

        The minimum feedback edge set of a digraph is a set of edges that
        intersect all the circuits of the digraph.  Equivalently, a minimum
        feedback arc set of a DiGraph is a set `S` of arcs such that the digraph
        `G-S` is acyclic. For more information, see the `Wikipedia article on
        feedback arc sets <http://en.wikipedia.org/wiki/Feedback_arc_set>`_.

        INPUT:

        - ``value_only`` -- boolean (default: ``False``)

          - When set to ``True``, only the minimum cardinal of a minimum edge
            set is returned.

          - When set to ``False``, the ``Set`` of edges of a minimal edge set is
            returned.

        - ``constraint_generation`` (boolean) -- whether to use constraint
          generation when solving the Mixed Integer Linear Program (default:
          ``True``).

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        ALGORITHM:

        This problem is solved using Linear Programming, in two different
        ways. The first one is to solve the following formulation:

        .. MATH::

            \mbox{Minimize : }&\sum_{(u,v)\in G} b_{(u,v)}\\
            \mbox{Such that : }&\\
            &\forall (u,v)\in G, d_u-d_v+ n \cdot b_{(u,v)}\geq 0\\
            &\forall u\in G, 0\leq d_u\leq |G|\\

        An explanation:

        An acyclic digraph can be seen as a poset, and every poset has a linear
        extension. This means that in any acyclic digraph the vertices can be
        ordered with a total order `<` in such a way that if `(u,v)\in G`, then
        `u<v`.

        Thus, this linear program is built in order to assign to each vertex `v`
        a number `d_v\in [0,\dots,n-1]` such that if there exists an edge
        `(u,v)\in G` such that `d_v<d_u`, then the edge `(u,v)` is removed.

        The number of edges removed is then minimized, which is the objective.

        (Constraint Generation)

        If the parameter ``constraint_generation`` is enabled, a more efficient
        formulation is used :

        .. MATH::

            \mbox{Minimize : }&\sum_{(u,v)\in G} b_{(u,v)}\\
            \mbox{Such that : }&\\
            &\forall C\text{ circuits }\subseteq G, \sum_{uv\in C}b_{(u,v)}\geq 1\\

        As the number of circuits contained in a graph is exponential, this LP
        is solved through constraint generation. This means that the solver is
        sequentially asked to solved the problem, knowing only a portion of the
        circuits contained in `G`, each time adding to the list of its
        constraints the circuit which its last answer had left intact.

        EXAMPLES:

        If the digraph is created from a graph, and hence is symmetric (if `uv`
        is an edge, then `vu` is an edge too), then obviously the cardinality of
        its feedback arc set is the number of edges in the first graph::

            sage: cycle=graphs.CycleGraph(5)
            sage: dcycle=DiGraph(cycle)
            sage: cycle.size()
            5
            sage: dcycle.feedback_edge_set(value_only=True)
            5

        And in this situation, for any edge `uv` of the first graph, `uv` of
        `vu` is in the returned feedback arc set::

           sage: g = graphs.RandomGNP(5,.3)
           sage: dg = DiGraph(g)
           sage: feedback = dg.feedback_edge_set()
           sage: (u,v,l) = g.edge_iterator().next()
           sage: (u,v) in feedback or (v,u) in feedback
           True

        TESTS:

        Comparing with/without constraint generation. Also double-checks ticket :trac:`12833`::

            sage: for i in range(20):
            ...      g = digraphs.RandomDirectedGNP(10,.3)
            ...      x = g.feedback_edge_set(value_only = True)
            ...      y = g.feedback_edge_set(value_only = True,
            ...             constraint_generation = False)
            ...      if x != y:
            ...         print "Oh my, oh my !"
            ...         break
        """
        # It would be a pity to start a LP if the digraph is already acyclic
        if self.is_directed_acyclic():
            return 0 if value_only else []

        from sage.numerical.mip import MixedIntegerLinearProgram

        ########################################
        # Constraint Generation Implementation #
        ########################################
        if constraint_generation:

            p = MixedIntegerLinearProgram(constraint_generation = True,
                                          maximization = False)

            # An variable for each edge
            b = p.new_variable(binary = True)

            # Variables are binary, and their coefficient in the objective is 1

            p.set_objective( p.sum( b[u,v]
                                  for u,v in self.edges(labels = False)))

            p.solve(log = verbose)

            # For as long as we do not break because the digraph is
            # acyclic....
            while True:

                # Building the graph without the edges removed by the LP
                h = DiGraph()
                for u,v in self.edges(labels = False):
                    if p.get_values(b[u,v]) < .5:
                        h.add_edge(u,v)

                # Is the digraph acyclic ?
                isok, certificate = h.is_directed_acyclic(certificate = True)

                # If so, we are done !
                if isok:
                    break

                if verbose:
                    print "Adding a constraint on circuit : ",certificate

                # There is a circuit left. Let's add the corresponding
                # constraint !

                p.add_constraint(
                    p.sum( b[u,v] for u,v in
                         zip(certificate, certificate[1:] + [certificate[0]])),
                    min = 1)

                obj = p.solve(log = verbose)

            if value_only:
                return Integer(round(obj))

            else:

                # listing the edges contained in the MFAS
                return [(u,v) for u,v in self.edges(labels = False)
                        if p.get_values(b[u,v]) > .5]

        ######################################
        # Ordering-based MILP Implementation #
        ######################################
        else:
            p=MixedIntegerLinearProgram(maximization=False, solver=solver)

            b=p.new_variable(binary = True)
            d=p.new_variable(integer = True)

            n=self.order()

            for (u,v) in self.edges(labels=None):
                p.add_constraint(d[u]-d[v]+n*(b[(u,v)]),min=1)

            for v in self:
                p.add_constraint(d[v] <= n)

            p.set_objective(p.sum([b[(u,v)] for (u,v) in self.edges(labels=None)]))

            if value_only:
                return Integer(round(p.solve(objective_only=True, log=verbose)))
            else:
                p.solve(log=verbose)

                b_sol=p.get_values(b)

                return [(u,v) for (u,v) in self.edges(labels=None) if b_sol[(u,v)]==1]

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

    def reverse_edge(self, u, v=None, label=None, inplace=True, multiedges=None):
        """
        Reverses the edge from u to v.

        INPUT:

        - ``inplace`` -- (default: True) if ``False``, a new digraph is created
           and returned as output, otherwise ``self`` is modified.

        - ``multiedges`` -- (default: None) how to decide what should be done in
          case of doubt (for instance when edge `(1,2)` is to be reversed in a
          graph while `(2,1)` already exists).

           - If set to ``True``, input graph will be forced to allow parallel
             edges if necessary and edge `(1,2)` will appear twice in the graph.

           - If set to ``False``, only one edge `(1,2)` will remain in the graph
             after `(2,1)` is reversed. Besides, the label of edge `(1,2)` will
             be overwritten with the label of edge `(2,1)`.

           The default behaviour (``multiedges = None``) will raise an exception
           each time a subjective decision (setting ``multiedges`` to ``True``
           or ``False``) is necessary to perform the operation.

        The following forms are all accepted:

        - D.reverse_edge( 1, 2 )
        - D.reverse_edge( (1, 2) )
        - D.reverse_edge( [1, 2] )
        - D.reverse_edge( 1, 2, 'label' )
        - D.reverse_edge( ( 1, 2, 'label') )
        - D.reverse_edge( [1, 2, 'label'] )
        - D.reverse_edge( ( 1, 2), label='label') )

        EXAMPLES:

        If ``inplace`` is ``True`` (default value), ``self`` is modified::

            sage: D = DiGraph([(0,1,2)])
            sage: D.reverse_edge(0,1)
            sage: D.edges()
            [(1, 0, 2)]

        If ``inplace`` is ``False``, ``self`` is not modified
        and a new digraph is returned::

            sage: D = DiGraph([(0,1,2)])
            sage: re = D.reverse_edge(0,1, inplace=False)
            sage: re.edges()
            [(1, 0, 2)]
            sage: D.edges()
            [(0, 1, 2)]

        If ``multiedges`` is ``True``, ``self`` will be forced to allow parallel
        edges when and only when it is necessary::

            sage: D = DiGraph( [(1, 2, 'A'), (2, 1, 'A'), (2, 3, None)] )
            sage: D.reverse_edge(1,2, multiedges=True)
            sage: D.edges()
            [(2, 1, 'A'), (2, 1, 'A'), (2, 3, None)]
            sage: D.allows_multiple_edges()
            True

        Even if ``multiedges`` is ``True``, ``self`` will not be forced to allow
        parallel edges when it is not necessary::

            sage: D = DiGraph( [(1,2,'A'), (2,1,'A'), (2, 3, None)] )
            sage: D.reverse_edge(2,3, multiedges=True)
            sage: D.edges()
            [(1, 2, 'A'), (2, 1, 'A'), (3, 2, None)]
            sage: D.allows_multiple_edges()
            False

        If user specifies ``multiedges = False``, ``self`` will not be forced to
        allow parallel edges and a parallel edge will get deleted::

            sage: D = DiGraph( [(1, 2, 'A'), (2, 1,'A'), (2, 3, None)] )
            sage: D.edges()
            [(1, 2, 'A'), (2, 1, 'A'), (2, 3, None)]
            sage: D.reverse_edge(1,2, multiedges=False)
            sage: D.edges()
            [(2, 1, 'A'), (2, 3, None)]

        Note that in the following graph, specifying ``multiedges = False`` will
        result in overwriting the label of `(1,2)` with the label of `(2,1)`::

            sage: D = DiGraph( [(1, 2, 'B'), (2, 1,'A'), (2, 3, None)] )
            sage: D.edges()
            [(1, 2, 'B'), (2, 1, 'A'), (2, 3, None)]
            sage: D.reverse_edge(2,1, multiedges=False)
            sage: D.edges()
            [(1, 2, 'A'), (2, 3, None)]

        If input edge in digraph has weight/label, then the weight/label should
        be preserved in the output digraph.  User does not need to specify the
        weight/label when calling function::

            sage: D = DiGraph([[0,1,2],[1,2,1]], weighted=True)
            sage: D.reverse_edge(0,1)
            sage: D.edges()
            [(1, 0, 2), (1, 2, 1)]
            sage: re = D.reverse_edge([1,2],inplace=False)
            sage: re.edges()
            [(1, 0, 2), (2, 1, 1)]

        If ``self`` has multiple copies (parallel edges) of the input edge, only
        1 of the parallel edges is reversed::

            sage: D = DiGraph([(0,1,'01'),(0,1,'01'),(0,1,'cat'),(1,2,'12')], weighted = True, multiedges = true)
            sage: re = D.reverse_edge([0,1,'01'],inplace=False)
            sage: re.edges()
            [(0, 1, '01'), (0, 1, 'cat'), (1, 0, '01'), (1, 2, '12')]

        If ``self`` has multiple copies (parallel edges) of the input edge but
        with distinct labels and no input label is specified, only 1 of the
        parallel edges is reversed (the edge that is labeled by the first label
        on the list returned by :meth:`.edge_label`)::

            sage: D = DiGraph([(0,1,'A'),(0,1,'B'),(0,1,'mouse'),(0,1,'cat')], multiedges = true)
            sage: D.edge_label(0,1)
            ['cat', 'mouse', 'B', 'A']
            sage: D.reverse_edge(0,1)
            sage: D.edges()
            [(0, 1, 'A'), (0, 1, 'B'), (0, 1, 'mouse'), (1, 0, 'cat')]

        Finally, an exception is raised when Sage does not know how to chose
        between allowing multiple edges and losing some data::

            sage: D = DiGraph([(0,1,'A'),(1,0,'B')])
            sage: D.reverse_edge(0,1)
            Traceback (most recent call last):
            ...
            ValueError: Reversing the given edge is about to create two parallel
            edges but input digraph doesn't allow them - User needs to specify
            multiedges is True or False.

        The following syntax is supported, but note that you must use
        the ``label`` keyword::

            sage: D = DiGraph()
            sage: D.add_edge((1,2), label='label')
            sage: D.edges()
            [(1, 2, 'label')]
            sage: D.reverse_edge((1,2),label ='label')
            sage: D.edges()
            [(2, 1, 'label')]
            sage: D.add_edge((1,2),'label')
            sage: D.edges()
            [(2, 1, 'label'), ((1, 2), 'label', None)]
            sage: D.reverse_edge((1,2), 'label')
            sage: D.edges()
            [(2, 1, 'label'), ('label', (1, 2), None)]

        TESTS::

            sage: D = DiGraph([(0,1,None)])
            sage: D.reverse_edge(0,1,'mylabel')
            Traceback (most recent call last):
            ...
            ValueError: Input edge must exist in the digraph.
        """
        # Assigns the expected values to u,v, and label depending on the input.
        if label is None:
            if v is None:
                try:
                    u, v, label = u
                except Exception:
                    try:
                        u, v = u
                    except Exception:
                        pass
        else:
            if v is None:
                try:
                    u, v = u
                except Exception:
                    pass

        if not self.has_edge(u,v,label):
            raise ValueError("Input edge must exist in the digraph.")

        tempG = self if inplace else self.copy()

        if label is None:
            if not tempG.allows_multiple_edges():
                label = tempG.edge_label(u,v)
            else:
                # If digraph has parallel edges for input edge, pick the first
                # from the labels on the list
                label = tempG.edge_label(u,v)[0]

        if ((not tempG.allows_multiple_edges()) and (tempG.has_edge(v,u))):
            # If user wants to force digraph to allow parallel edges
            if multiedges == True:
                tempG.allow_multiple_edges(True)
                tempG.delete_edge(u,v,label)
                tempG.add_edge(v,u,label)

            # If user does not want to force digraph to allow parallel
            # edges, we delete edge u to v and overwrite v,u with the
            # label of u,v
            elif multiedges == False:
                tempG.delete_edge(u,v,label)
                tempG.set_edge_label(v,u,label)

            # User is supposed to specify multiedges True or None
            else:
                raise ValueError("Reversing the given edge is about to "
                                 "create two parallel edges but input digraph "
                                 "doesn't allow them - User needs to specify "
                                 "multiedges is True or False.")
        else:
            tempG.delete_edge(u,v,label)
            tempG.add_edge(v,u,label)

        if not inplace:
            return tempG

    def reverse_edges(self, edges, inplace=True, multiedges=None):
        """
        Reverses a list of edges.

        INPUT:

        - ``edges`` -- a list of edges in the DiGraph.

        - ``inplace`` -- (default: True) if ``False``, a new digraph is created
           and returned as output, otherwise ``self`` is modified.

        - ``multiedges`` -- (default: None) if ``True``, input graph will be
           forced to allow parallel edges when necessary (for more information
           see the documentation of :meth:`~DiGraph.reverse_edge`)

        .. SEEALSO::

            :meth:`~DiGraph.reverse_edge` - Reverses a single edge.

        EXAMPLES:

        If ``inplace`` is ``True`` (default value), ``self`` is modified::

            sage: D = DiGraph({ 0: [1,1,3], 2: [3,3], 4: [1,5]}, multiedges = true)
            sage: D.reverse_edges( [ [0,1], [0,3] ])
            sage: D.reverse_edges( [ (2,3),(4,5) ])
            sage: D.edges()
            [(0, 1, None), (1, 0, None), (2, 3, None), (3, 0, None),
             (3, 2, None), (4, 1, None), (5, 4, None)]

        If ``inplace`` is ``False``, ``self`` is not modified and a new digraph
        is returned::

            sage: D = DiGraph ([(0,1,'A'),(1,0,'B'),(1,2,'C')])
            sage: re = D.reverse_edges( [ (0,1), (1,2) ],
            ...                         inplace = False,
            ...                         multiedges = True)
            sage: re.edges()
            [(1, 0, 'A'), (1, 0, 'B'), (2, 1, 'C')]
            sage: D.edges()
            [(0, 1, 'A'), (1, 0, 'B'), (1, 2, 'C')]
            sage: D.allows_multiple_edges()
            False
            sage: re.allows_multiple_edges()
            True

        If ``multiedges`` is ``True``, ``self`` will be forced to allow parallel
        edges when and only when it is necessary::

            sage: D = DiGraph( [(1, 2, 'A'), (2, 1, 'A'), (2, 3, None)] )
            sage: D.reverse_edges([(1,2),(2,3)], multiedges=True)
            sage: D.edges()
            [(2, 1, 'A'), (2, 1, 'A'), (3, 2, None)]
            sage: D.allows_multiple_edges()
            True

        Even if ``multiedges`` is ``True``, ``self`` will not be forced to allow
        parallel edges when it is not necessary::

            sage: D = DiGraph( [(1, 2, 'A'), (2, 1, 'A'), (2,3, None)] )
            sage: D.reverse_edges([(2,3)], multiedges=True)
            sage: D.edges()
            [(1, 2, 'A'), (2, 1, 'A'), (3, 2, None)]
            sage: D.allows_multiple_edges()
            False

        If ``multiedges`` is ``False``, ``self`` will not be forced to allow
        parallel edges and an edge will get deleted::

            sage: D = DiGraph( [(1,2), (2,1)] )
            sage: D.edges()
            [(1, 2, None), (2, 1, None)]
            sage: D.reverse_edges([(1,2)], multiedges=False)
            sage: D.edges()
            [(2, 1, None)]

        If input edge in digraph has weight/label, then the weight/label should
        be preserved in the output digraph.  User does not need to specify the
        weight/label when calling function::

            sage: D = DiGraph([(0,1,'01'),(1,2,1),(2,3,'23')], weighted = True)
            sage: D.reverse_edges([(0,1,'01'),(1,2),(2,3)])
            sage: D.edges()
            [(1, 0, '01'), (2, 1, 1), (3, 2, '23')]

        TESTS::

            sage: D = digraphs.Circuit(6)
            sage: D.reverse_edges(D.edges(),inplace=False).edges()
            [(0, 5, None), (1, 0, None), (2, 1, None),
             (3, 2, None), (4, 3, None), (5, 4, None)]

            sage: D = digraphs.Kautz(2,3)
            sage: Dr = D.reverse_edges(D.edges(),inplace=False,multiedges=True)
            sage: Dr.edges() == D.reverse().edges()
            True
        """
        tempG = self if inplace else self.copy()
        for e in edges:
            tempG.reverse_edge(e,inplace=True,multiedges=multiedges)
        if not inplace:
            return tempG

    ### Paths and cycles iterators

    def _all_paths_iterator(self, vertex, ending_vertices=None,
                            simple=False, max_length=None, trivial=False):
        r"""
        Returns an iterator over the paths of self starting with the
        given vertex.

        INPUT:

        -  ``vertex`` - the starting vertex of the paths.
        -  ``ending_vertices`` - iterable (default: None) on the allowed
           ending vertices of the paths. If None, then all vertices are
           allowed.
        -  ``simple`` - boolean (default: False). If set to True, then
           only simple paths are considered. Simple paths are paths in
           which no two arcs share a head or share a tail, i.e. every
           vertex in the path is entered at most once and exited at most
           once.
        -  ``max_length`` - non negative integer (default: None). The
           maximum length of the enumerated paths. If set to None, then
           all lengths are allowed.
        -  ``trivial`` - boolean (default: False). If set to True, then
           the empty paths are also enumerated.

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

        One may wish to enumerate simple paths, which are paths in which
        no two arcs share a head or share a tail, i.e. every vertex in
        the path is entered at most once and exited at most once. The
        result is always finite but may take a long time to compute::

            sage: pi = g._all_paths_iterator('a', simple=True)
            sage: list(pi)
            [['a', 'a'], ['a', 'b'], ['a', 'b', 'c'], ['a', 'b', 'c', 'd']]
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
            [['a'], ['a', 'a'], ['a', 'b'], ['a', 'a', 'a'], ['a', 'a', 'b'],
             ['a', 'b', 'c'], ['a', 'a', 'a', 'a'], ['a', 'a', 'a', 'b'],
             ['a', 'a', 'b', 'c'], ['a', 'b', 'c', 'd']]
        """
        if ending_vertices is None:
            ending_vertices = self
        if max_length is None:
            from sage.rings.infinity import Infinity
            max_length = Infinity
        if max_length < 1:
            return

        # Start with the empty path; we will try all extensions of it
        queue = []
        path = [vertex]

        if trivial and vertex in ending_vertices:
            yield path
        while True:
            # Build next generation of paths, one arc longer; max_length refers
            # to edges and not vertices, hence <= and not <
            if len(path) <= max_length:

                # We try all possible extensions
                if simple:
                    # We only keep simple extensions. An extension is simple
                    # iff the new vertex being entered has not previously
                    # occurred in the path, or has occurred but only been
                    # exited (i.e. is the first vertex in the path). In this
                    # latter case we must not exit the new vertex again, so we
                    # do not consider it for further extension, but just yield
                    # it immediately. See trac #12385.
                    for neighbor in self.neighbor_out_iterator(path[-1]):
                        if neighbor not in path:
                            queue.append(path + [neighbor])
                        elif ( neighbor == path[0] and
                               neighbor in ending_vertices ):
                            yield path + [neighbor]

                else:
                    # Non-simple paths requested: we add all of them
                    for neighbor in self.neighbor_out_iterator(path[-1]):
                        queue.append(path + [neighbor])

            if not queue:
                break
            path = queue.pop(0)     # get the next path

            if path[-1] in ending_vertices:
                yield path      # yield good path


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
           then only simple paths are considered. These are paths in
           which no two arcs share a head or share a tail, i.e. every
           vertex in the path is entered at most once and exited at most
           once.
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
        :meth:`all_simple_paths`)::

            sage: pi = g.all_paths_iterator(simple=True)
            sage: list(pi)
            [['a', 'a'], ['a', 'b'], ['b', 'c'], ['c', 'd'], ['d', 'c'],
             ['a', 'b', 'c'], ['b', 'c', 'd'], ['c', 'd', 'c'],
             ['d', 'c', 'd'], ['a', 'b', 'c', 'd']]

        Or simply bound the length of the enumerated paths::

            sage: pi = g.all_paths_iterator(starting_vertices=['a'], ending_vertices=['b', 'c'], max_length=6)
            sage: list(pi)
            [['a', 'b'], ['a', 'a', 'b'], ['a', 'b', 'c'],
             ['a', 'a', 'a', 'b'], ['a', 'a', 'b', 'c'],
             ['a', 'a', 'a', 'a', 'b'], ['a', 'a', 'a', 'b', 'c'],
             ['a', 'b', 'c', 'd', 'c'], ['a', 'a', 'a', 'a', 'a', 'b'],
             ['a', 'a', 'a', 'a', 'b', 'c'], ['a', 'a', 'b', 'c', 'd', 'c'],
             ['a', 'a', 'a', 'a', 'a', 'a', 'b'],
             ['a', 'a', 'a', 'a', 'a', 'b', 'c'],
             ['a', 'a', 'a', 'b', 'c', 'd', 'c'],
             ['a', 'b', 'c', 'd', 'c', 'd', 'c']]

        By default, empty paths are not enumerated, but it may be
        parametrized::

            sage: pi = g.all_paths_iterator(simple=True, trivial=True)
            sage: list(pi)
            [['a'], ['b'], ['c'], ['d'], ['a', 'a'], ['a', 'b'], ['b', 'c'],
             ['c', 'd'], ['d', 'c'], ['a', 'b', 'c'], ['b', 'c', 'd'],
             ['c', 'd', 'c'], ['d', 'c', 'd'], ['a', 'b', 'c', 'd']]
            sage: pi = g.all_paths_iterator(simple=True, trivial=False)
            sage: list(pi)
            [['a', 'a'], ['a', 'b'], ['b', 'c'], ['c', 'd'], ['d', 'c'],
             ['a', 'b', 'c'], ['b', 'c', 'd'], ['c', 'd', 'c'],
             ['d', 'c', 'd'], ['a', 'b', 'c', 'd']]
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
        with one of the given vertices. Simple paths are paths in which
        no two arcs share a head or share a tail, i.e. every vertex in
        the path is entered at most once and exited at most once.

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
            [['a', 'a'], ['a', 'b'], ['b', 'c'], ['c', 'd'], ['d', 'c'],
             ['a', 'b', 'c'], ['b', 'c', 'd'], ['c', 'd', 'c'],
             ['d', 'c', 'd'], ['a', 'b', 'c', 'd']]

        One may compute all paths having specific starting and/or
        ending vertices::

            sage: g.all_simple_paths(starting_vertices=['a'])
            [['a', 'a'], ['a', 'b'], ['a', 'b', 'c'], ['a', 'b', 'c', 'd']]
            sage: g.all_simple_paths(starting_vertices=['a'], ending_vertices=['c'])
            [['a', 'b', 'c']]
            sage: g.all_simple_paths(starting_vertices=['a'], ending_vertices=['b', 'c'])
            [['a', 'b'], ['a', 'b', 'c']]

        It is also possible to bound the length of the paths::

            sage: g.all_simple_paths(max_length=2)
            [['a', 'a'], ['a', 'b'], ['b', 'c'], ['c', 'd'], ['d', 'c'],
             ['a', 'b', 'c'], ['b', 'c', 'd'], ['c', 'd', 'c'],
             ['d', 'c', 'd']]

        By default, empty paths are not enumerated, but this can
        be parametrized::

            sage: g.all_simple_paths(starting_vertices=['a'], trivial=True)
            [['a'], ['a', 'a'], ['a', 'b'], ['a', 'b', 'c'],
             ['a', 'b', 'c', 'd']]
            sage: g.all_simple_paths(starting_vertices=['a'], trivial=False)
            [['a', 'a'], ['a', 'b'], ['a', 'b', 'c'], ['a', 'b', 'c', 'd']]
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

            See also :meth:`all_simple_cycles`.

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
            [['a', 'a'], ['a', 'a', 'a'], ['c', 'd', 'c'],
             ['a', 'a', 'a', 'a']]

        By default, cycles differing only by their starting point are not all
        enumerated, but this may be parametrized::

            sage: it = g.all_cycles_iterator(max_length=3, rooted=False)
            sage: list(it)
            [['a', 'a'], ['a', 'a', 'a'], ['c', 'd', 'c'],
             ['a', 'a', 'a', 'a']]
            sage: it = g.all_cycles_iterator(max_length=3, rooted=True)
            sage: list(it)
            [['a', 'a'], ['a', 'a', 'a'], ['c', 'd', 'c'], ['d', 'c', 'd'],
             ['a', 'a', 'a', 'a']]

        One may prefer to enumerate simple cycles, i.e. cycles such that the only
        vertex occuring twice in it is the starting and ending one (see also
        :meth:`all_simple_cycles`)::

            sage: it = g.all_cycles_iterator(simple=True)
            sage: list(it)
            [['a', 'a'], ['c', 'd', 'c']]
            sage: g = digraphs.Circuit(4)
            sage: list(g.all_cycles_iterator(simple=True))
            [[0, 1, 2, 3, 0]]
        """
        if starting_vertices is None:
            starting_vertices = self
        # Since a cycle is always included in a given strongly connected
        # component, we may remove edges from the graph
        sccs = self.strongly_connected_components()
        d = {}
        for id, component in enumerate(sccs):
            for v in component:
                d[v] = id
        h = self.copy()
        h.delete_edges([ (u,v) for (u,v) in h.edge_iterator(labels=False)
                if d[u] != d[v] ])
        # We create one cycles iterator per vertex. This is necessary if we
        # want to iterate over cycles with increasing length.
        vertex_iterators = dict([(v, h._all_cycles_iterator_vertex( v
                                        , starting_vertices=starting_vertices
                                        , simple=simple
                                        , rooted=rooted
                                        , max_length=max_length
                                        , trivial=trivial
                                        , remove_acyclic_edges=False
                                        )) for v in starting_vertices])
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
            # We update the cycle iterator to its next available cycle if it
            # exists
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

            Although the number of simple cycles of a finite graph is
            always finite, computing all its cycles may take a very long
            time.

        EXAMPLES::

            sage: g = DiGraph({'a' : ['a', 'b'], 'b' : ['c'], 'c' : ['d'], 'd' : ['c']}, loops=True)
            sage: g.all_simple_cycles()
            [['a', 'a'], ['c', 'd', 'c']]

        The directed version of the Petersen graph::

            sage: g = graphs.PetersenGraph().to_directed()
            sage: g.all_simple_cycles(max_length=4)
            [[0, 1, 0], [0, 4, 0], [0, 5, 0], [1, 2, 1], [1, 6, 1], [2, 3, 2],
             [2, 7, 2], [3, 8, 3], [3, 4, 3], [4, 9, 4], [5, 8, 5], [5, 7, 5],
             [6, 8, 6], [6, 9, 6], [7, 9, 7]]
            sage: g.all_simple_cycles(max_length=6)
            [[0, 1, 0], [0, 4, 0], [0, 5, 0], [1, 2, 1], [1, 6, 1], [2, 3, 2],
             [2, 7, 2], [3, 8, 3], [3, 4, 3], [4, 9, 4], [5, 8, 5], [5, 7, 5],
             [6, 8, 6], [6, 9, 6], [7, 9, 7], [0, 1, 2, 3, 4, 0],
             [0, 1, 2, 7, 5, 0], [0, 1, 6, 8, 5, 0], [0, 1, 6, 9, 4, 0],
             [0, 4, 9, 6, 1, 0], [0, 4, 9, 7, 5, 0], [0, 4, 3, 8, 5, 0],
             [0, 4, 3, 2, 1, 0], [0, 5, 8, 3, 4, 0], [0, 5, 8, 6, 1, 0],
             [0, 5, 7, 9, 4, 0], [0, 5, 7, 2, 1, 0], [1, 2, 3, 8, 6, 1],
             [1, 2, 7, 9, 6, 1], [1, 6, 8, 3, 2, 1], [1, 6, 9, 7, 2, 1],
             [2, 3, 8, 5, 7, 2], [2, 3, 4, 9, 7, 2], [2, 7, 9, 4, 3, 2],
             [2, 7, 5, 8, 3, 2], [3, 8, 6, 9, 4, 3], [3, 4, 9, 6, 8, 3],
             [5, 8, 6, 9, 7, 5], [5, 7, 9, 6, 8, 5], [0, 1, 2, 3, 8, 5, 0],
             [0, 1, 2, 7, 9, 4, 0], [0, 1, 6, 8, 3, 4, 0],
             [0, 1, 6, 9, 7, 5, 0], [0, 4, 9, 6, 8, 5, 0],
             [0, 4, 9, 7, 2, 1, 0], [0, 4, 3, 8, 6, 1, 0],
             [0, 4, 3, 2, 7, 5, 0], [0, 5, 8, 3, 2, 1, 0],
             [0, 5, 8, 6, 9, 4, 0], [0, 5, 7, 9, 6, 1, 0],
             [0, 5, 7, 2, 3, 4, 0], [1, 2, 3, 4, 9, 6, 1],
             [1, 2, 7, 5, 8, 6, 1], [1, 6, 8, 5, 7, 2, 1],
             [1, 6, 9, 4, 3, 2, 1], [2, 3, 8, 6, 9, 7, 2],
             [2, 7, 9, 6, 8, 3, 2], [3, 8, 5, 7, 9, 4, 3],
             [3, 4, 9, 7, 5, 8, 3]]

        The complete graph (without loops) on `4` vertices::

            sage: g = graphs.CompleteGraph(4).to_directed()
            sage: g.all_simple_cycles()
            [[0, 1, 0], [0, 2, 0], [0, 3, 0], [1, 2, 1], [1, 3, 1], [2, 3, 2],
             [0, 1, 2, 0], [0, 1, 3, 0], [0, 2, 1, 0], [0, 2, 3, 0],
             [0, 3, 1, 0], [0, 3, 2, 0], [1, 2, 3, 1], [1, 3, 2, 1],
             [0, 1, 2, 3, 0], [0, 1, 3, 2, 0], [0, 2, 1, 3, 0],
             [0, 2, 3, 1, 0], [0, 3, 1, 2, 0], [0, 3, 2, 1, 0]]

        If the graph contains a large number of cycles, one can bound
        the length of the cycles, or simply restrict the possible
        starting vertices of the cycles::

            sage: g = graphs.CompleteGraph(20).to_directed()
            sage: g.all_simple_cycles(max_length=2)
            [[0, 1, 0], [0, 2, 0], [0, 3, 0], [0, 4, 0], [0, 5, 0], [0, 6, 0],
             [0, 7, 0], [0, 8, 0], [0, 9, 0], [0, 10, 0], [0, 11, 0],
             [0, 12, 0], [0, 13, 0], [0, 14, 0], [0, 15, 0], [0, 16, 0],
             [0, 17, 0], [0, 18, 0], [0, 19, 0], [1, 2, 1], [1, 3, 1],
             [1, 4, 1], [1, 5, 1], [1, 6, 1], [1, 7, 1], [1, 8, 1], [1, 9, 1],
             [1, 10, 1], [1, 11, 1], [1, 12, 1], [1, 13, 1], [1, 14, 1],
             [1, 15, 1], [1, 16, 1], [1, 17, 1], [1, 18, 1], [1, 19, 1],
             [2, 3, 2], [2, 4, 2], [2, 5, 2], [2, 6, 2], [2, 7, 2], [2, 8, 2],
             [2, 9, 2], [2, 10, 2], [2, 11, 2], [2, 12, 2], [2, 13, 2],
             [2, 14, 2], [2, 15, 2], [2, 16, 2], [2, 17, 2], [2, 18, 2],
             [2, 19, 2], [3, 4, 3], [3, 5, 3], [3, 6, 3], [3, 7, 3], [3, 8, 3],
             [3, 9, 3], [3, 10, 3], [3, 11, 3], [3, 12, 3], [3, 13, 3],
             [3, 14, 3], [3, 15, 3], [3, 16, 3], [3, 17, 3], [3, 18, 3],
             [3, 19, 3], [4, 5, 4], [4, 6, 4], [4, 7, 4], [4, 8, 4], [4, 9, 4],
             [4, 10, 4], [4, 11, 4], [4, 12, 4], [4, 13, 4], [4, 14, 4],
             [4, 15, 4], [4, 16, 4], [4, 17, 4], [4, 18, 4], [4, 19, 4],
             [5, 6, 5], [5, 7, 5], [5, 8, 5], [5, 9, 5], [5, 10, 5],
             [5, 11, 5], [5, 12, 5], [5, 13, 5], [5, 14, 5], [5, 15, 5],
             [5, 16, 5], [5, 17, 5], [5, 18, 5], [5, 19, 5], [6, 7, 6],
             [6, 8, 6], [6, 9, 6], [6, 10, 6], [6, 11, 6], [6, 12, 6],
             [6, 13, 6], [6, 14, 6], [6, 15, 6], [6, 16, 6], [6, 17, 6],
             [6, 18, 6], [6, 19, 6], [7, 8, 7], [7, 9, 7], [7, 10, 7],
             [7, 11, 7], [7, 12, 7], [7, 13, 7], [7, 14, 7], [7, 15, 7],
             [7, 16, 7], [7, 17, 7], [7, 18, 7], [7, 19, 7], [8, 9, 8],
             [8, 10, 8], [8, 11, 8], [8, 12, 8], [8, 13, 8], [8, 14, 8],
             [8, 15, 8], [8, 16, 8], [8, 17, 8], [8, 18, 8], [8, 19, 8],
             [9, 10, 9], [9, 11, 9], [9, 12, 9], [9, 13, 9], [9, 14, 9],
             [9, 15, 9], [9, 16, 9], [9, 17, 9], [9, 18, 9], [9, 19, 9],
             [10, 11, 10], [10, 12, 10], [10, 13, 10], [10, 14, 10],
             [10, 15, 10], [10, 16, 10], [10, 17, 10], [10, 18, 10],
             [10, 19, 10], [11, 12, 11], [11, 13, 11], [11, 14, 11],
             [11, 15, 11], [11, 16, 11], [11, 17, 11], [11, 18, 11],
             [11, 19, 11], [12, 13, 12], [12, 14, 12], [12, 15, 12],
             [12, 16, 12], [12, 17, 12], [12, 18, 12], [12, 19, 12],
             [13, 14, 13], [13, 15, 13], [13, 16, 13], [13, 17, 13],
             [13, 18, 13], [13, 19, 13], [14, 15, 14], [14, 16, 14],
             [14, 17, 14], [14, 18, 14], [14, 19, 14], [15, 16, 15],
             [15, 17, 15], [15, 18, 15], [15, 19, 15], [16, 17, 16],
             [16, 18, 16], [16, 19, 16], [17, 18, 17], [17, 19, 17],
             [18, 19, 18]]
            sage: g = graphs.CompleteGraph(20).to_directed()
            sage: g.all_simple_cycles(max_length=2, starting_vertices=[0])
            [[0, 1, 0], [0, 2, 0], [0, 3, 0], [0, 4, 0], [0, 5, 0], [0, 6, 0],
             [0, 7, 0], [0, 8, 0], [0, 9, 0], [0, 10, 0], [0, 11, 0],
             [0, 12, 0], [0, 13, 0], [0, 14, 0], [0, 15, 0], [0, 16, 0],
             [0, 17, 0], [0, 18, 0], [0, 19, 0]]

        One may prefer to distinguish equivalent cycles having distinct
        starting vertices (compare the following examples)::

            sage: g = graphs.CompleteGraph(4).to_directed()
            sage: g.all_simple_cycles(max_length=2, rooted=False)
            [[0, 1, 0], [0, 2, 0], [0, 3, 0], [1, 2, 1], [1, 3, 1], [2, 3, 2]]
            sage: g.all_simple_cycles(max_length=2, rooted=True)
            [[0, 1, 0], [0, 2, 0], [0, 3, 0], [1, 0, 1], [1, 2, 1], [1, 3, 1],
             [2, 0, 2], [2, 1, 2], [2, 3, 2], [3, 0, 3], [3, 1, 3], [3, 2, 3]]
        """
        return list(self.all_cycles_iterator(starting_vertices=starting_vertices, simple=True, rooted=rooted, max_length=max_length, trivial=trivial))

    def path_semigroup(self):
        """
        The partial semigroup formed by the paths of this quiver.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','c']}, 2:{3:['b']}})
            sage: F = Q.path_semigroup(); F
            Partial semigroup formed by the directed paths of Multi-digraph on 3 vertices
            sage: list(F)
            [e_1, e_2, e_3, a, c, b, a*b, c*b]

        """
        from sage.quivers.path_semigroup import PathSemigroup
        # If self is immutable, then the copy is really cheap:
        # __copy__ just returns self.
        return PathSemigroup(self.copy(immutable=True))

    ### Directed Acyclic Graphs (DAGs)

    def topological_sort(self, implementation = "default"):
        """
        Returns a topological sort of the digraph if it is acyclic, and
        raises a TypeError if the digraph contains a directed cycle. As
        topological sorts are not necessarily unique, different
        implementations may yield different results.

        A topological sort is an ordering of the vertices of the digraph
        such that each vertex comes before all of its successors. That
        is, if `u` comes before `v` in the sort, then there may be
        a directed path from `u` to `v`, but there will be no directed
        path from `v` to `u`.

        INPUT:

        - ``implementation`` -- Use the default Cython implementation
          (``implementation = default``), the default NetworkX library
          (``implementation = "NetworkX"``) or the recursive NetworkX
          implementation (``implementation = "recursive"``)

        .. SEEALSO::

            - :meth:`is_directed_acyclic` -- Tests whether a directed
              graph is acyclic (can also join a certificate --
              a topological sort or a circuit in the graph1).

        EXAMPLES::

            sage: D = DiGraph({ 0:[1,2,3], 4:[2,5], 1:[8], 2:[7], 3:[7],
            ...     5:[6,7], 7:[8], 6:[9], 8:[10], 9:[10] })
            sage: D.plot(layout='circular').show()
            sage: D.topological_sort()
            [4, 5, 6, 9, 0, 1, 2, 3, 7, 8, 10]

        ::

            sage: D.add_edge(9,7)
            sage: D.topological_sort()
            [4, 5, 6, 9, 0, 1, 2, 3, 7, 8, 10]

        Using the NetworkX implementation ::

            sage: D.topological_sort(implementation = "NetworkX")
            [4, 5, 6, 9, 0, 1, 2, 3, 7, 8, 10]

        Using the NetworkX recursive implementation ::

            sage: D.topological_sort(implementation = "recursive")
            [4, 5, 6, 9, 0, 3, 2, 7, 1, 8, 10]

        ::

            sage: D.add_edge(7,4)
            sage: D.topological_sort()
            Traceback (most recent call last):
            ...
            TypeError: Digraph is not acyclic-- there is no topological
            sort.

        .. note::

           There is a recursive version of this in NetworkX, it used to
           have problems in earlier versions but they have since been
           fixed::

              sage: import networkx
              sage: D = DiGraph({ 0:[1,2,3], 4:[2,5], 1:[8], 2:[7], 3:[7],
              ...     5:[6,7], 7:[8], 6:[9], 8:[10], 9:[10] })
              sage: N = D.networkx_graph()
              sage: networkx.topological_sort(N)
              [4, 5, 6, 9, 0, 1, 2, 3, 7, 8, 10]
              sage: networkx.topological_sort_recursive(N)
              [4, 5, 6, 9, 0, 3, 2, 7, 1, 8, 10]

        TESTS:

        A wrong value for the ``implementation`` keyword::

            sage: D.topological_sort(implementation = "cloud-reading")
            Traceback (most recent call last):
            ...
            ValueError: implementation must be set to one of "default"
            or "NetworkX"
        """

        if implementation == "default":
            b, ordering = self._backend.is_directed_acyclic(certificate = True)
            if b:
                return ordering
            else:
                raise TypeError('Digraph is not acyclic-- there is no topological sort.')

        elif implementation == "NetworkX" or implementation == "recursive":
            import networkx
            if implementation == "NetworkX":
                S = networkx.topological_sort(self.networkx_graph(copy=False))
            else:
                S = networkx.topological_sort_recursive(self.networkx_graph(copy=False))
            if S is None:
                raise TypeError('Digraph is not acyclic-- there is no topological sort.')
            else:
                return S

        else:
            raise ValueError("implementation must be set to one of \"default\" or \"NetworkX\"")

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

    def layout_acyclic(self, **options):
        """
        Computes a ranked layout so that all edges point upward.

        To this end, the heights of the vertices are set according to the level
        set decomposition of the graph (see :meth:`.level_sets`).

        This is achieved by calling ``graphviz`` and ``dot2tex`` if available
        (see :meth:`.layout_graphviz`), and using a random horizontal
        placement of the vertices otherwise (see :meth:`.layout_acyclic_dummy`).

        Non acyclic graphs are partially supported by ``graphviz``, which then
        chooses some edges to point down.

        EXAMPLES::

            sage: H = DiGraph({0:[1,2],1:[3],2:[3],3:[],5:[1,6],6:[2,3]})

        The actual layout computed will depend on whether dot2tex and
        graphviz are installed, so we don't test it.

            sage: H.layout_acyclic()
            {0: [..., ...], 1: [..., ...], 2: [..., ...], 3: [..., ...], 5: [..., ...], 6: [..., ...]}

        """
        if have_dot2tex():
            return self.layout_graphviz(**options)
        else:
            return self.layout_acyclic_dummy(**options)

    def layout_acyclic_dummy(self, heights = None, **options):
        """
        Computes a (dummy) ranked layout of an acyclic graph so that
        all edges point upward. To this end, the heights of the
        vertices are set according to the level set decomposition of
        the graph (see :meth:`.level_sets`).

        EXAMPLES::

            sage: H = DiGraph({0:[1,2],1:[3],2:[3],3:[],5:[1,6],6:[2,3]})
            sage: H.layout_acyclic_dummy()
            {0: [1.00..., 0], 1: [1.00..., 1], 2: [1.51..., 2], 3: [1.50..., 3], 5: [2.01..., 0], 6: [2.00..., 1]}

            sage: H = DiGraph({0:[1,2],1:[3],2:[3],3:[1],5:[1,6],6:[2,3]})
            sage: H.layout_acyclic_dummy()
            Traceback (most recent call last):
            ...
            AssertionError: `self` should be an acyclic graph
        """
        if heights is None:
            assert self.is_directed_acyclic(), "`self` should be an acyclic graph"
            levels = self.level_sets()
            levels = [sorted(z) for z in levels]
            heights = dict([[i, levels[i]] for i in range(len(levels))])
        return self.layout_ranked(heights = heights, **options)

    def level_sets(self):
        """
        Returns the level set decomposition of the digraph.

        OUTPUT:

         - a list of non empty lists of vertices of this graph

        The level set decomposition of the digraph is a list `l` such that the
        level `l[i]` contains all the vertices having all their predecessors in
        the levels `l[j]` for `j<i`, and at least one in level `l[i-1]` (unless
        `i=0`).

        The level decomposition contains exactly the vertices not occuring in
        any cycle of the graph. In particular, the graph is acyclic if and only
        if the decomposition forms a set partition of its vertices, and we
        recover the usual level set decomposition of the corresponding poset.

        EXAMPLES::

            sage: H = DiGraph({0:[1,2],1:[3],2:[3],3:[],5:[1,6],6:[2,3]})
            sage: H.level_sets()
            [[0, 5], [1, 6], [2], [3]]

            sage: H = DiGraph({0:[1,2],1:[3],2:[3],3:[1],5:[1,6],6:[2,3]})
            sage: H.level_sets()
            [[0, 5], [6], [2]]

        This routine is mostly used for Hasse diagrams of posets::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[]})
            sage: [len(x) for x in H.level_sets()]
            [1, 2, 1]

        ::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2], 1:[3], 2:[4], 3:[4]})
            sage: [len(x) for x in H.level_sets()]
            [1, 2, 1, 1]

        Complexity: `O(n+m)` in time and `O(n)` in memory (besides the
        storage of the graph itself), where `n` and `m` are
        respectively the number of vertices and edges (assuming that
        appending to a list is constant time, which it is not quite).
        """
        in_degrees = self.in_degree(labels=True)
        level = [x for x in in_degrees if in_degrees[x]==0]
        Levels = []
        while len(level) != 0:
            Levels.append(level)
            new_level = []
            for x in level:
                for y in self.neighbors_out(x):
                    in_degrees[y] -= 1
                    if in_degrees[y] == 0:
                        new_level.append(y)
            level = new_level
        return Levels

    def strongly_connected_components(self):
        """
        Returns the list of strongly connected components.

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

        TESTS:

        Checking against NetworkX, and another of Sage's implementations::

            sage: from sage.graphs.base.static_sparse_graph import strongly_connected_components
            sage: import networkx
            sage: for i in range(100):                                     # long
            ...        g = digraphs.RandomDirectedGNP(100,.05)             # long
            ...        h = g.networkx_graph()                              # long
            ...        scc1 = g.strongly_connected_components()            # long
            ...        scc2 = networkx.strongly_connected_components(h)    # long
            ...        scc3 = strongly_connected_components(g)             # long
            ...        s1 = Set(map(Set,scc1))                             # long
            ...        s2 = Set(map(Set,scc2))                             # long
            ...        s3 = Set(map(Set,scc3))                             # long
            ...        if s1 != s2:                                        # long
            ...            print "Ooch !"                                  # long
            ...        if s1 != s3:                                        # long
            ...            print "Oooooch !"                               # long

        """

        try:
            vertices = set(self.vertices())
            scc = []
            while vertices:
                tmp = self.strongly_connected_component_containing_vertex(vertices.__iter__().next())
                vertices.difference_update(set(tmp))
                scc.append(tmp)
            return scc

        except AttributeError:
            import networkx
            return networkx.strongly_connected_components(self.networkx_graph(copy=False))

    def strongly_connected_component_containing_vertex(self, v):
        """
        Returns the strongly connected component containing a given vertex

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
            raise AttributeError("This function is only defined for C graphs.")

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

    def strongly_connected_components_digraph(self, keep_labels = False):
        r"""
        Returns the digraph of the strongly connected components

        INPUT:

         - ``keep_labels`` -- boolean (default: False)

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

        The following digraph has three strongly connected components,
        and the digraph of those is a chain::

            sage: g = DiGraph({0:{1:"01", 2: "02", 3: 03}, 1: {2: "12"}, 2:{1: "21", 3: "23"}})
            sage: scc_digraph = g.strongly_connected_components_digraph()
            sage: scc_digraph.vertices()
            [{0}, {3}, {1, 2}]
            sage: scc_digraph.edges()
            [({0}, {3}, None), ({0}, {1, 2}, None), ({1, 2}, {3}, None)]

        By default, the labels are discarded, and the result has no
        loops nor multiple edges. If ``keep_labels`` is ``True``, then
        the labels are kept, and the result is a multi digraph,
        possibly with multiple edges and loops. However, edges in the
        result with same source, target, and label are not duplicated
        (see the edges from 0 to the strongly connected component
        `\{1,2\}` below)::

            sage: g = DiGraph({0:{1:"0-12", 2: "0-12", 3: "0-3"}, 1: {2: "1-2", 3: "1-3"}, 2:{1: "2-1", 3: "2-3"}})
            sage: scc_digraph = g.strongly_connected_components_digraph(keep_labels = True)
            sage: scc_digraph.vertices()
            [{0}, {3}, {1, 2}]
            sage: scc_digraph.edges()
            [({0}, {3}, '0-3'), ({0}, {1, 2}, '0-12'),
             ({1, 2}, {3}, '1-3'), ({1, 2}, {3}, '2-3'),
             ({1, 2}, {1, 2}, '1-2'), ({1, 2}, {1, 2}, '2-1')]

        """

        from sage.sets.set import Set

        scc = self.strongly_connected_components()
        scc_set = map(Set,scc)

        d = {}
        for i,c in enumerate(scc):
            for v in c:
                d[v] = i

        if keep_labels:
            g = DiGraph(multiedges=True, loops=True)
            g.add_vertices(range(len(scc)))

            g.add_edges( set((d[u], d[v], label) for (u,v,label) in self.edges() ) )
            g.relabel(scc_set, inplace = True)

        else:
            g = DiGraph(multiedges=False, loops=False)
            g.add_vertices(range(len(scc)))

            for u,v in self.edges(labels=False):
                g.add_edge(d[u], d[v])

            g.relabel(scc_set, inplace = True)

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

    def is_aperiodic(self):
        r"""
        Return whether the current ``DiGraph`` is aperiodic.

        A directed graph is aperiodic if there is no integer ``k > 1``
        that divides the length of every cycle in the graph, cf.
        :wikipedia:`Aperiodic_graph`.

        EXAMPLES:

        The following graph has period ``2``, so it is not aperiodic::

            sage: g = DiGraph({ 0: [1], 1: [0] })
            sage: g.is_aperiodic()
            False

        The following graph has a cycle of length 2 and a cycle of length 3,
        so it is aperiodic::

            sage: g = DiGraph({ 0: [1, 4], 1: [2], 2: [0], 4: [0]})
            sage: g.is_aperiodic()
            True

        """
        import networkx
        return networkx.is_aperiodic(self.networkx_graph(copy=False))

import types

import sage.graphs.comparability
DiGraph.is_transitive = types.MethodType(sage.graphs.comparability.is_transitive, None, DiGraph)
