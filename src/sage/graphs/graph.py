r"""
Graph Theory

This module implements many graph theoretic operations and
concepts.

AUTHORS:

-  Robert L. Miller (2006-10-22): initial version

-  William Stein (2006-12-05): Editing

-  Robert L. Miller (2007-01-13): refactoring, adjusting for
   NetworkX-0.33, fixed plotting bugs (2007-01-23): basic tutorial,
   edge labels, loops, multiple edges and arcs (2007-02-07): graph6
   and sparse6 formats, matrix input

-  Emily Kirkmann (2007-02-11): added graph_border option to plot
   and show

-  Robert L. Miller (2007-02-12): vertex color-maps, graph
   boundaries, graph6 helper functions in Cython

-  Robert L. Miller Sage Days 3 (2007-02-17-21): 3d plotting in
   Tachyon

-  Robert L. Miller (2007-02-25): display a partition

-  Robert L. Miller (2007-02-28): associate arbitrary objects to
   vertices, edge and arc label display (in 2d), edge coloring

-  Robert L. Miller (2007-03-21): Automorphism group, isomorphism
   check, canonical label

-  Robert L. Miller (2007-06-07-09): NetworkX function wrapping

-  Michael W. Hansen (2007-06-09): Topological sort generation

-  Emily Kirkman, Robert L. Miller Sage Days 4: Finished wrapping
   NetworkX

-  Emily Kirkman (2007-07-21): Genus (including circular planar,
   all embeddings and all planar embeddings), all paths, interior
   paths

-  Bobby Moretti (2007-08-12): fixed up plotting of graphs with
   edge colors differentiated by label

-  Jason Grout (2007-09-25): Added functions, bug fixes, and
   general enhancements

-  Robert L. Miller (Sage Days 7): Edge labeled graph isomorphism

-  Tom Boothby (Sage Days 7): Miscellaneous awesomeness

-  Tom Boothby (2008-01-09): Added graphviz output

-  David Joyner (2009-2): Fixed docstring bug related to GAP.

-  Stephen Hartke (2009-07-26): Fixed bug in blocks_and_cut_vertices()
   that caused an incorrect result when the vertex 0 was a cut vertex.

-  Stephen Hartke (2009-08-22): Fixed bug in blocks_and_cut_vertices()
   where the list of cut_vertices is not treated as a set.

-  Anders Jonsson (2009-10-10): Counting of spanning trees and out-trees added.

-  Nathann Cohen (2009-09) : Cliquer, Connectivity, Flows
                             and everything that uses Linear Programming
                             and class numerical.MIP


Graph Format
------------

The Sage Graph Class: NetworkX plus
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sage graphs are actually NetworkX graphs, wrapped in a Sage class.
In fact, any graph can produce its underlying NetworkX graph. For
example,

::

    sage: import networkx
    sage: G = graphs.PetersenGraph()
    sage: N = G.networkx_graph()
    sage: isinstance(N, networkx.graph.Graph)
    True

The NetworkX graph is essentially a dictionary of dictionaries::

    sage: N.adj
    {0: {1: None, 4: None, 5: None}, 1: {0: None, 2: None, 6: None}, 2: {1: None, 3: None, 7: None}, 3: {8: None, 2: None, 4: None}, 4: {0: None, 9: None, 3: None}, 5: {0: None, 8: None, 7: None}, 6: {8: None, 1: None, 9: None}, 7: {9: None, 2: None, 5: None}, 8: {3: None, 5: None, 6: None}, 9: {4: None, 6: None, 7: None}}

Each dictionary key is a vertex label, and each key in the
following dictionary is a neighbor of that vertex. In undirected
graphs, there is redundancy: for example, the dictionary containing
the entry ``1: {2: None}`` implies it must contain
``{2: {1: None}``. The innermost entry of ``None`` is
related to edge labeling (see section :ref:`Graph:labels`).

Supported formats
~~~~~~~~~~~~~~~~~

Sage Graphs can be created from a wide range of inputs. A few
examples are covered here.


-  NetworkX dictionary format:

   ::

       sage: d = {0: [1,4,5], 1: [2,6], 2: [3,7], 3: [4,8], 4: [9], \
             5: [7, 8], 6: [8,9], 7: [9]}
       sage: G = Graph(d); G
       Graph on 10 vertices
       sage: G.plot().show()    # or G.show()

-  A NetworkX graph:

   ::

       sage: K = networkx.complete_bipartite_graph(12,7)
       sage: G = Graph(K)
       sage: G.degree()
       [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 12, 12, 12, 12, 12, 12, 12]

-  graph6 or sparse6 format:

   ::

       sage: s = ':I`AKGsaOs`cI]Gb~'
       sage: G = Graph(s, sparse=True); G
       Looped multi-graph on 10 vertices
       sage: G.plot().show()    # or G.show()

   Note that the ``\`` character is an escape character in Python, and
   also a character used by graph6 strings:

   ::

       sage: G = Graph('Ihe\n@GUA')
       Traceback (most recent call last):
       ...
       RuntimeError: The string (Ihe) seems corrupt: for n = 10, the string is too short.

   In Python, the escaped character ``\`` is represented by ``\\``:

   ::

       sage: G = Graph('Ihe\\n@GUA')
       sage: G.plot().show()    # or G.show()

-  adjacency matrix: In an adjacency matrix, each column and each
   row represent a vertex. If a 1 shows up in row `i`, column
   `j`, there is an edge `(i,j)`.

   ::

       sage: M = Matrix([(0,1,0,0,1,1,0,0,0,0),(1,0,1,0,0,0,1,0,0,0), \
       (0,1,0,1,0,0,0,1,0,0), (0,0,1,0,1,0,0,0,1,0),(1,0,0,1,0,0,0,0,0,1), \
       (1,0,0,0,0,0,0,1,1,0), (0,1,0,0,0,0,0,0,1,1),(0,0,1,0,0,1,0,0,0,1), \
       (0,0,0,1,0,1,1,0,0,0), (0,0,0,0,1,0,1,1,0,0)])
       sage: M
       [0 1 0 0 1 1 0 0 0 0]
       [1 0 1 0 0 0 1 0 0 0]
       [0 1 0 1 0 0 0 1 0 0]
       [0 0 1 0 1 0 0 0 1 0]
       [1 0 0 1 0 0 0 0 0 1]
       [1 0 0 0 0 0 0 1 1 0]
       [0 1 0 0 0 0 0 0 1 1]
       [0 0 1 0 0 1 0 0 0 1]
       [0 0 0 1 0 1 1 0 0 0]
       [0 0 0 0 1 0 1 1 0 0]
       sage: G = Graph(M); G
       Graph on 10 vertices
       sage: G.plot().show()    # or G.show()

-  incidence matrix: In an incidence matrix, each row represents a
   vertex and each column represents an edge.

   ::

       sage: M = Matrix([(-1,0,0,0,1,0,0,0,0,0,-1,0,0,0,0), \
       (1,-1,0,0,0,0,0,0,0,0,0,-1,0,0,0),(0,1,-1,0,0,0,0,0,0,0,0,0,-1,0,0), \
       (0,0,1,-1,0,0,0,0,0,0,0,0,0,-1,0),(0,0,0,1,-1,0,0,0,0,0,0,0,0,0,-1), \
       (0,0,0,0,0,-1,0,0,0,1,1,0,0,0,0),(0,0,0,0,0,0,0,1,-1,0,0,1,0,0,0), \
       (0,0,0,0,0,1,-1,0,0,0,0,0,1,0,0),(0,0,0,0,0,0,0,0,1,-1,0,0,0,1,0), \
       (0,0,0,0,0,0,1,-1,0,0,0,0,0,0,1)])
       sage: M
       [-1  0  0  0  1  0  0  0  0  0 -1  0  0  0  0]
       [ 1 -1  0  0  0  0  0  0  0  0  0 -1  0  0  0]
       [ 0  1 -1  0  0  0  0  0  0  0  0  0 -1  0  0]
       [ 0  0  1 -1  0  0  0  0  0  0  0  0  0 -1  0]
       [ 0  0  0  1 -1  0  0  0  0  0  0  0  0  0 -1]
       [ 0  0  0  0  0 -1  0  0  0  1  1  0  0  0  0]
       [ 0  0  0  0  0  0  0  1 -1  0  0  1  0  0  0]
       [ 0  0  0  0  0  1 -1  0  0  0  0  0  1  0  0]
       [ 0  0  0  0  0  0  0  0  1 -1  0  0  0  1  0]
       [ 0  0  0  0  0  0  1 -1  0  0  0  0  0  0  1]
       sage: G = Graph(M); G
       Graph on 10 vertices
       sage: G.plot().show()    # or G.show()
       sage: DiGraph(matrix(2,[0,0,-1,1]), format="incidence_matrix")
       Traceback (most recent call last):
       ...
       ValueError: There must be two nonzero entries (-1 & 1) per column.


Generators
----------

If you wish to iterate through all the isomorphism types of graphs,
type, for example::

    sage: for g in graphs(4):
    ...     print g.spectrum()
    [0, 0, 0, 0]
    [1, 0, 0, -1]
    [1.4142135623..., 0, 0, -1.4142135623...]
    [2, 0, -1, -1]
    [1.7320508075..., 0, 0, -1.7320508075...]
    [1, 1, -1, -1]
    [1.6180339887..., 0.6180339887..., -0.6180339887..., -1.6180339887...]
    [2.1700864866..., 0.3111078174..., -1, -1.4811943040...]
    [2, 0, 0, -2]
    [2.5615528128..., 0, -1, -1.5615528128...]
    [3, -1, -1, -1]

For some commonly used graphs to play with, type

::

    sage: graphs.[tab]          # not tested

and hit {tab}. Most of these graphs come with their own custom
plot, so you can see how people usually visualize these graphs.

::

    sage: G = graphs.PetersenGraph()
    sage: G.plot().show()    # or G.show()
    sage: G.degree_histogram()
    [0, 0, 0, 10]
    sage: G.adjacency_matrix()
    [0 1 0 0 1 1 0 0 0 0]
    [1 0 1 0 0 0 1 0 0 0]
    [0 1 0 1 0 0 0 1 0 0]
    [0 0 1 0 1 0 0 0 1 0]
    [1 0 0 1 0 0 0 0 0 1]
    [1 0 0 0 0 0 0 1 1 0]
    [0 1 0 0 0 0 0 0 1 1]
    [0 0 1 0 0 1 0 0 0 1]
    [0 0 0 1 0 1 1 0 0 0]
    [0 0 0 0 1 0 1 1 0 0]

::

    sage: S = G.subgraph([0,1,2,3])
    sage: S.plot().show()    # or S.show()
    sage: S.density()
    1/2

::

    sage: G = GraphQuery(display_cols=['graph6'], num_vertices=7, diameter=5)
    sage: L = G.get_graphs_list()
    sage: graphs_list.show_graphs(L)

.. _Graph:labels:

Labels
------

Each vertex can have any hashable object as a label. These are
things like strings, numbers, and tuples. Each edge is given a
default label of ``None``, but if specified, edges can
have any label at all. Edges between vertices `u` and
`v` are represented typically as ``(u, v, l)``, where
``l`` is the label for the edge.

Note that vertex labels themselves cannot be mutable items::

    sage: M = Matrix( [[0,0],[0,0]] )
    sage: G = Graph({ 0 : { M : None } })
    Traceback (most recent call last):
    ...
    TypeError: mutable matrices are unhashable

However, if one wants to define a dictionary, with the same keys
and arbitrary objects for entries, one can make that association::

    sage: d = {0 : graphs.DodecahedralGraph(), 1 : graphs.FlowerSnark(), \
          2 : graphs.MoebiusKantorGraph(), 3 : graphs.PetersenGraph() }
    sage: d[2]
    Moebius-Kantor Graph: Graph on 16 vertices
    sage: T = graphs.TetrahedralGraph()
    sage: T.vertices()
    [0, 1, 2, 3]
    sage: T.set_vertices(d)
    sage: T.get_vertex(1)
    Flower Snark: Graph on 20 vertices

Database
--------

There is a database available for searching for graphs that satisfy
a certain set of parameters, including number of vertices and
edges, density, maximum and minimum degree, diameter, radius, and
connectivity. To see a list of all search parameter keywords broken
down by their designated table names, type

::

    sage: graph_db_info()
    {...}

For more details on data types or keyword input, enter

::

    sage: GraphQuery?    # not tested

The results of a query can be viewed with the show method, or can be
viewed individually by iterating through the results:

::

    sage: Q = GraphQuery(display_cols=['graph6'],num_vertices=7, diameter=5)
    sage: Q.show()
    Graph6
    --------------------
    F@?]O
    F@OKg
    F?`po
    F?gqg
    FIAHo
    F@R@o
    FA_pW
    FGC{o
    FEOhW

Show each graph as you iterate through the results:

::

    sage: for g in Q:
    ...     show(g)


Visualization
-------------

To see a graph `G` you are working with, there
are three main options. You can view the graph in two dimensions via
matplotlib with ``show()``. ::

    sage: G = graphs.RandomGNP(15,.3)
    sage: G.show()

And you can view it in three dimensions via jmol with ``show3d()``. ::

    sage: G.show3d()

Or it can be rendered with `\mbox{\rm\LaTeX}`.  This requires the right
additions to a standard `\mbox{\rm\TeX}` installation.  Then standard
Sage commands, such as ``view(G)`` will display the graph, or
``latex(G)`` will produce a string suitable for inclusion in a
`\mbox{\rm\LaTeX}` document.  More details on this are at
the :mod:`sage.graphs.graph_latex` module. ::

    sage: from sage.graphs.graph_latex import check_tkz_graph
    sage: check_tkz_graph()  # random - depends on TeX installation
    sage: latex(G)
    \begin{tikzpicture}
    ...
    \end{tikzpicture}

Graph classes and methods
-------------------------
"""

#*****************************************************************************
#      Copyright (C) 2006 - 2007 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer import Integer

import sage.graphs.generic_graph_pyx as generic_graph_pyx
from sage.graphs.generic_graph import GenericGraph
from sage.graphs.digraph import DiGraph

class Graph(GenericGraph):
    r"""
    Undirected graph.

    A graph is a set of vertices connected by edges
    (cf. http://en.wikipedia.org/wiki/Graph_(mathematics) ).

    One can very easily create a graph in Sage by typing::

        sage: g = Graph()

    By typing the name of the graph, one can get some basic information
    about it::

        sage: g
        Graph on 0 vertices

    This graph is not very interesting as it is by default the empty graph.
    But Sage contains a large collection of pre-defined graph classes that
    can be listed this way:

    * Within a Sage session, type ``graphs.``
      (Do not press "Enter", and do not forget the final period ".")
    * Hit "tab".

    You will see a list of methods which will construct named graphs. For
    example::

        sage: g = graphs.PetersenGraph()
        sage: g.plot()

    or::

        sage: g = graphs.ChvatalGraph()
        sage: g.plot()

    In order to obtain more information about these graph constructors, access
    the documentation using the command ``graphs.RandomGNP?``.

    Once you have defined the graph you want, you can begin to work on it
    by using the almost 200 functions on graphs in the Sage library!
    If your graph is named ``g``, you can list these functions as previously
    this way

    * Within a Sage session, type ``g.``
      (Do not press "Enter", and do not forget the final period "." )
    * Hit "tab".

    As usual, you can get some information about what these functions do by
    typing (e.g. if you want to know about the ``diameter()`` method)
    ``g.diameter?``.

    If you have defined a graph ``g`` having several connected components
    (i.e. ``g.is_connected()`` returns False), you can print each one of its
    connected components with only two lines::

        sage: for component in g.connected_components():
        ...      g.subgraph(component).plot()


    INPUT:

    -  ``data`` -- can be any of the following:

      #.  A dictionary of dictionaries

      #.  A dictionary of lists

      #.  A NumPy matrix or ndarray

      #.  A Sage adjacency matrix or incidence matrix

      #.  A pygraphviz agraph

      #.  A SciPy sparse matrix

      #.  A NetworkX digraph

    -  ``pos`` -  a positioning dictionary: for example, the
       spring layout from NetworkX for the 5-cycle is::

         {0: [-0.91679746, 0.88169588],
          1: [ 0.47294849, 1.125     ],
          2: [ 1.125     ,-0.12867615],
          3: [ 0.12743933,-1.125     ],
          4: [-1.125     ,-0.50118505]}

    -  ``name`` - (must be an explicitly named parameter,
       i.e., ``name="complete")`` gives the graph a name

    -  ``loops`` - boolean, whether to allow loops (ignored
       if data is an instance of the ``Graph`` class)

    -  ``multiedges`` - boolean, whether to allow multiple
       edges (ignored if data is an instance of the ``Graph`` class)

    -  ``weighted`` - whether graph thinks of itself as
       weighted or not. See ``self.weighted()``

    -  ``format`` - if None, Graph tries to guess- can be
       several values, including:

       -  ``'graph6'`` - Brendan McKay's graph6 format, in a
          string (if the string has multiple graphs, the first graph is
          taken)

       -  ``'sparse6'`` - Brendan McKay's sparse6 format, in a
          string (if the string has multiple graphs, the first graph is
          taken)

       -  ``'adjacency_matrix'`` - a square Sage matrix M,
          with M[i,j] equal to the number of edges {i,j}

       -  ``'weighted_adjacency_matrix'`` - a square Sage
          matrix M, with M[i,j] equal to the weight of the single edge {i,j}.
          Given this format, weighted is ignored (assumed True).

       -  ``'incidence_matrix'`` - a Sage matrix, with one
          column C for each edge, where if C represents {i, j}, C[i] is -1
          and C[j] is 1

       -  ``'elliptic_curve_congruence'`` - data must be an
          iterable container of elliptic curves, and the graph produced has
          each curve as a vertex (it's Cremona label) and an edge E-F
          labelled p if and only if E is congruent to F mod p

    -  ``boundary`` - a list of boundary vertices, if
       empty, graph is considered as a 'graph without boundary'

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

    We illustrate the first six input formats (the other two
    involve packages that are currently not standard in Sage):

    #. A dictionary of dictionaries::

        sage: g = Graph({0:{1:'x',2:'z',3:'a'}, 2:{5:'out'}}); g
        Graph on 5 vertices

       The labels ('x', 'z', 'a', 'out') are labels for edges. For
       example, 'out' is the label for the edge on 2 and 5. Labels can be
       used as weights, if all the labels share some common parent.

       ::

        sage: a,b,c,d,e,f = sorted(SymmetricGroup(3))
        sage: Graph({b:{d:'c',e:'p'}, c:{d:'p',e:'c'}})
        Graph on 4 vertices

    #. A dictionary of lists::

        sage: g = Graph({0:[1,2,3], 2:[4]}); g
        Graph on 5 vertices

    #. A list of vertices and a function describing adjacencies. Note
       that the list of vertices and the function must be enclosed in a
       list (i.e., [list of vertices, function]).

       Construct the Paley graph over GF(13).

       ::

          sage: g=Graph([GF(13), lambda i,j: i!=j and (i-j).is_square()])
          sage: g.vertices()
          [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
          sage: g.adjacency_matrix()
          [0 1 0 1 1 0 0 0 0 1 1 0 1]
          [1 0 1 0 1 1 0 0 0 0 1 1 0]
          [0 1 0 1 0 1 1 0 0 0 0 1 1]
          [1 0 1 0 1 0 1 1 0 0 0 0 1]
          [1 1 0 1 0 1 0 1 1 0 0 0 0]
          [0 1 1 0 1 0 1 0 1 1 0 0 0]
          [0 0 1 1 0 1 0 1 0 1 1 0 0]
          [0 0 0 1 1 0 1 0 1 0 1 1 0]
          [0 0 0 0 1 1 0 1 0 1 0 1 1]
          [1 0 0 0 0 1 1 0 1 0 1 0 1]
          [1 1 0 0 0 0 1 1 0 1 0 1 0]
          [0 1 1 0 0 0 0 1 1 0 1 0 1]
          [1 0 1 1 0 0 0 0 1 1 0 1 0]

       Construct the line graph of a complete graph.

       ::

          sage: g=graphs.CompleteGraph(4)
          sage: line_graph=Graph([g.edges(labels=false), \
                 lambda i,j: len(set(i).intersection(set(j)))>0], \
                 loops=False)
          sage: line_graph.vertices()
          [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
          sage: line_graph.adjacency_matrix()
          [0 1 1 1 1 0]
          [1 0 1 1 0 1]
          [1 1 0 0 1 1]
          [1 1 0 0 1 1]
          [1 0 1 1 0 1]
          [0 1 1 1 1 0]

    #. A NumPy matrix or ndarray::

        sage: import numpy
        sage: A = numpy.array([[0,1,1],[1,0,1],[1,1,0]])
        sage: Graph(A)
        Graph on 3 vertices

    #. A graph6 or sparse6 string: Sage automatically recognizes
       whether a string is in graph6 or sparse6 format::

           sage: s = ':I`AKGsaOs`cI]Gb~'
           sage: Graph(s,sparse=True)
           Looped multi-graph on 10 vertices

       ::

           sage: G = Graph('G?????')
           sage: G = Graph("G'?G?C")
           Traceback (most recent call last):
           ...
           RuntimeError: The string seems corrupt: valid characters are
           ?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
           sage: G = Graph('G??????')
           Traceback (most recent call last):
           ...
           RuntimeError: The string (G??????) seems corrupt: for n = 8, the string is too long.

       ::

          sage: G = Graph(":I'AKGsaOs`cI]Gb~")
          Traceback (most recent call last):
          ...
          RuntimeError: The string seems corrupt: valid characters are
          ?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~

       There are also list functions to take care of lists of graphs::

           sage: s = ':IgMoqoCUOqeb\n:I`AKGsaOs`cI]Gb~\n:I`EDOAEQ?PccSsge\N\n'
           sage: graphs_list.from_sparse6(s)
           [Looped multi-graph on 10 vertices, Looped multi-graph on 10 vertices, Looped multi-graph on 10 vertices]

    #. A Sage matrix:
       Note: If format is not specified, then Sage assumes a symmetric square
       matrix is an adjacency matrix, otherwise an incidence matrix.

       - an adjacency matrix::

            sage: M = graphs.PetersenGraph().am(); M
            [0 1 0 0 1 1 0 0 0 0]
            [1 0 1 0 0 0 1 0 0 0]
            [0 1 0 1 0 0 0 1 0 0]
            [0 0 1 0 1 0 0 0 1 0]
            [1 0 0 1 0 0 0 0 0 1]
            [1 0 0 0 0 0 0 1 1 0]
            [0 1 0 0 0 0 0 0 1 1]
            [0 0 1 0 0 1 0 0 0 1]
            [0 0 0 1 0 1 1 0 0 0]
            [0 0 0 0 1 0 1 1 0 0]
            sage: Graph(M)
            Graph on 10 vertices

         ::

            sage: Graph(matrix([[1,2],[2,4]]),loops=True,sparse=True)
            Looped multi-graph on 2 vertices

            sage: M = Matrix([[0,1,-1],[1,0,-1/2],[-1,-1/2,0]]); M
            [   0    1   -1]
            [   1    0 -1/2]
            [  -1 -1/2    0]
            sage: G = Graph(M,sparse=True); G
            Graph on 3 vertices
            sage: G.weighted()
            True

       - an incidence matrix::

            sage: M = Matrix(6, [-1,0,0,0,1, 1,-1,0,0,0, 0,1,-1,0,0, 0,0,1,-1,0, 0,0,0,1,-1, 0,0,0,0,0]); M
            [-1  0  0  0  1]
            [ 1 -1  0  0  0]
            [ 0  1 -1  0  0]
            [ 0  0  1 -1  0]
            [ 0  0  0  1 -1]
            [ 0  0  0  0  0]
            sage: Graph(M)
            Graph on 6 vertices

            sage: Graph(Matrix([[1],[1],[1]]))
            Traceback (most recent call last):
            ...
            ValueError: Non-symmetric or non-square matrix assumed to be an incidence matrix: There must be two nonzero entries (-1 & 1) per column.
            sage: Graph(Matrix([[1],[1],[0]]))
            Traceback (most recent call last):
            ...
            ValueError: Non-symmetric or non-square matrix assumed to be an incidence matrix: Each column represents an edge: -1 goes to 1.

            sage: M = Matrix([[0,1,-1],[1,0,-1],[-1,-1,0]]); M
            [ 0  1 -1]
            [ 1  0 -1]
            [-1 -1  0]
            sage: Graph(M,sparse=True)
            Graph on 3 vertices

            sage: M = Matrix([[0,1,1],[1,0,0],[0,0,0]]); M
            [0 1 1]
            [1 0 0]
            [0 0 0]
            sage: Graph(M)
            Traceback (most recent call last):
            ...
            ValueError: Non-symmetric or non-square matrix assumed to be an incidence matrix: There must be two nonzero entries (-1 & 1) per column.

            sage: M = Matrix([[0,1,1],[1,0,1],[-1,-1,0]]); M
            [ 0  1  1]
            [ 1  0  1]
            [-1 -1  0]
            sage: Graph(M)
            Traceback (most recent call last):
            ...
            ValueError: Non-symmetric or non-square matrix assumed to be an incidence matrix: Each column represents an edge: -1 goes to 1.


    #. A NetworkX XGraph::

          sage: import networkx
          sage: g = networkx.XGraph({0:[1,2,3], 2:[4]})
          sage: Graph(g)
          Graph on 5 vertices

    #. A NetworkX graph::

           sage: import networkx
           sage: g = networkx.Graph({0:[1,2,3], 2:[4]})
           sage: DiGraph(g)
           Digraph on 5 vertices

    Note that in all cases, we copy the NetworkX structure.

       ::

          sage: import networkx
          sage: g = networkx.Graph({0:[1,2,3], 2:[4]})
          sage: G = Graph(g, implementation='networkx')
          sage: H = Graph(g, implementation='networkx')
          sage: G._backend._nxg is H._backend._nxg
          False
    """
    _directed = False

    def __init__(self, data=None, pos=None, loops=None, format=None,
                 boundary=[], weighted=None, implementation='c_graph',
                 sparse=True, vertex_labels=True, **kwds):
        """
        TESTS::

            sage: G = Graph()
            sage: loads(dumps(G)) == G
            True
            sage: a = matrix(2,2,[1,0,0,1])
            sage: Graph(a).adjacency_matrix() == a
            True

            sage: a = matrix(2,2,[2,0,0,1])
            sage: Graph(a,sparse=True).adjacency_matrix() == a
            True
        """
        GenericGraph.__init__(self)
        msg = ''
        multiedges = kwds.get('multiedges', None)
        from sage.structure.element import is_Matrix
        from sage.misc.misc import uniq
        if format is None and isinstance(data, str):
            if data[:10] == ">>graph6<<":
                data = data[10:]
                format = 'graph6'
            elif data[:11] == ">>sparse6<<":
                data = data[11:]
                format = 'sparse6'
            elif data[0] == ':':
                format = 'sparse6'
            else:
                format = 'graph6'
        if format is None and is_Matrix(data):
            if data.is_square() and data == data.transpose():
                format = 'adjacency_matrix'
            else:
                format = 'incidence_matrix'
                msg += "Non-symmetric or non-square matrix assumed to be an incidence matrix: "
        if format is None and isinstance(data, Graph):
            format = 'Graph'
        from sage.graphs.all import DiGraph
        if format is None and isinstance(data, DiGraph):
            data = data.to_undirected()
            format = 'Graph'
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
            if isinstance(data, (networkx.DiGraph, networkx.XDiGraph)):
                data = data.to_undirected()
                format = 'NX'
            elif isinstance(data, (networkx.Graph, networkx.XGraph)):
                format = 'NX'
        if format is None and isinstance(data, (int, Integer)):
            format = 'int'
        if format is None and data is None:
            format = 'int'
            data = 0
        if format is None:
            import networkx
            data = networkx.XGraph(data)
            format = 'NX'

        # At this point, format has been set.

        verts = None

        if format == 'graph6':
            if loops      is None: loops      = False
            if weighted   is None: weighted   = False
            if multiedges is None: multiedges = False
            if not isinstance(data, str):
                raise ValueError('If input format is graph6, then data must be a string.')
            n = data.find('\n')
            if n == -1:
                n = len(data)
            ss = data[:n]
            n, s = generic_graph_pyx.N_inverse(ss)
            m = generic_graph_pyx.R_inverse(s, n)
            expected = n*(n-1)/2 + (6 - n*(n-1)/2)%6
            if len(m) > expected:
                raise RuntimeError("The string (%s) seems corrupt: for n = %d, the string is too long."%(ss,n))
            elif len(m) < expected:
                raise RuntimeError("The string (%s) seems corrupt: for n = %d, the string is too short."%(ss,n))
            num_verts = n
        elif format == 'sparse6':
            if loops      is None: loops      = True
            if weighted   is None: weighted   = False
            if multiedges is None: multiedges = True
            from math import ceil, floor
            from sage.misc.functional import log
            n = data.find('\n')
            if n == -1:
                n = len(data)
            s = data[:n]
            n, s = generic_graph_pyx.N_inverse(s[1:])
            if n == 0:
                edges = []
            else:
                k = int(ceil(log(n,2)))
                ords = [ord(i) for i in s]
                if any(o > 126 or o < 63 for o in ords):
                    raise RuntimeError("The string seems corrupt: valid characters are \n" + ''.join([chr(i) for i in xrange(63,127)]))
                bits = ''.join([generic_graph_pyx.binary(o-63).zfill(6) for o in ords])
                b = []
                x = []
                for i in xrange(int(floor(len(bits)/(k+1)))):
                    b.append(int(bits[(k+1)*i:(k+1)*i+1],2))
                    x.append(int(bits[(k+1)*i+1:(k+1)*i+k+1],2))
                v = 0
                edges = []
                for i in xrange(len(b)):
                    if b[i] == 1:
                        v += 1
                    if x[i] > v:
                        v = x[i]
                    else:
                        if v < n:
                            edges.append((x[i],v))
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
                        raise ValueError("Non-weighted graph's"+
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
                        raise ValueError("Non-looped graph's adjacency"+
                        " matrix must have zeroes on the diagonal.")
                    break
            num_verts = data.nrows()
        elif format == 'incidence_matrix':
            try:
                positions = []
                for c in data.columns():
                    NZ = c.nonzero_positions()
                    positions.append(tuple(NZ))
                    if len(NZ) != 2:
                        msg += "There must be two nonzero entries (-1 & 1) per column."
                        assert False
                    L = uniq(c.list())
                    L.sort()
                    if L != [-1,0,1]:
                        msg += "Each column represents an edge: -1 goes to 1."
                        assert False
                if loops      is None: loops     = False
                if weighted   is None: weighted  = False
                if multiedges is None:
                    total = len(positions)
                    multiedges = (  len(uniq(positions)) < total  )
            except AssertionError:
                raise ValueError(msg)
            num_verts = data.nrows()
        elif format == 'Graph':
            if loops is None: loops = data.allows_loops()
            elif not loops and data.has_loops():
                raise ValueError("No loops but input graph has loops.")
            if multiedges is None: multiedges = data.allows_multiple_edges()
            elif not multiedges:
                e = data.edges(labels=False)
                e = [sorted(f) for f in e]
                if len(e) != len(uniq(e)):
                    raise ValueError("No multiple edges but input graph"+
                    " has multiple edges.")
            if weighted is None: weighted = data.weighted()
            num_verts = data.num_verts()
            verts = data.vertex_iterator()
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
                    if hash(u) > hash(v):
                        if v in data and u in data[v]:
                            if data[u][v] != data[v][u]:
                                raise ValueError("Dict does not agree on edge (%s,%s)"%(u,v))
                            continue
                    if multiedges is not False and not isinstance(data[u][v], list):
                        if multiedges is None: multiedges = False
                        if multiedges:
                            raise ValueError("Dict of dicts for multigraph must be in the format {v : {u : list}}")
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
                verts=verts.union([v for v in data[u] if v not in verts])
                if len(uniq(data[u])) != len(data[u]):
                    if multiedges is False:
                        raise ValueError("Non-multigraph input dict has multiple edges (%s,%s)"%(u, choice([v for v in data[u] if data[u].count(v) > 1])))
                    if multiedges is None: multiedges = True
            if multiedges is None: multiedges = False
            num_verts = len(verts)
        elif format == 'NX':
            if weighted is None:
                if isinstance(data, networkx.Graph):
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
        elif format in ['int', 'elliptic_curve_congruence']:
            if weighted   is None: weighted   = False
            if multiedges is None: multiedges = False
            if loops      is None: loops      = False
            if format == 'int': num_verts = data
            else:
                num_verts = len(data)
                curves = data
                verts = [curve.cremona_label() for curve in data]

        # weighted, multiedges, loops, verts and num_verts should now be set

        if implementation == 'networkx':
            import networkx
            from sage.graphs.base.graph_backends import NetworkXGraphBackend
            if format == 'Graph':
                self._backend = NetworkXGraphBackend(data.networkx_graph())
                self._weighted = weighted
                self.allow_loops(loops)
                self.allow_multiple_edges(multiedges)
            else:
                self._backend = NetworkXGraphBackend(networkx.XGraph())
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
            if format == 'Graph':
                self._backend = CGB(0, directed=False)
                self.add_vertices(verts)
                self._weighted = weighted
                self.allow_loops(loops, check=False)
                self.allow_multiple_edges(multiedges, check=False)
                for u,v,l in data.edge_iterator():
                    self._backend.add_edge(u,v,l,False)
            else:
                if verts is not None:
                    self._backend = CGB(0, directed=False)
                    self.add_vertices(verts)
                else:
                    self._backend = CGB(num_verts, directed=False)
                self._weighted = weighted
                self.allow_loops(loops, check=False)
                self.allow_multiple_edges(multiedges, check=False)
            self._backend.directed = False
        else:
            raise NotImplementedError("Supported implementations: networkx, c_graph.")

        if format == 'graph6':
            k = 0
            for i in xrange(n):
                for j in xrange(i):
                    if m[k] == '1':
                        self.add_edge(i, j)
                    k += 1
        elif format == 'sparse6':
            for i,j in edges:
                self.add_edge(i,j)
        elif format == 'adjacency_matrix':
            e = []
            if weighted:
                for i,j in data.nonzero_positions():
                    if i <= j:
                        e.append((i,j,data[i][j]))
            elif multiedges:
                for i,j in data.nonzero_positions():
                    if i <= j:
                        e += [(i,j)]*int(data[i][j])
            else:
                for i,j in data.nonzero_positions():
                    if i <= j:
                        e.append((i,j))
            self.add_edges(e)
        elif format == 'incidence_matrix':
            self.add_edges(positions)
        elif format == 'Graph':
            self.name(data.name())
        elif format == 'rule':
            verts = list(verts)
            for u in xrange(num_verts):
                for v in xrange(u+1):
                    uu,vv = verts[u], verts[v]
                    if f(uu,vv):
                        self.add_edge(uu,vv)
        elif format == 'dict_of_dicts':
            for u in data:
                for v in data[u]:
                    if hash(u) <= hash(v) or v not in data or u not in data[v]:
                        if multiedges:
                            self.add_edges([(u,v,l) for l in data[u][v]])
                        else:
                            self.add_edge((u,v,data[u][v]))
        elif format == 'dict_of_lists':
            for u in data:
                for v in data[u]:
                    if multiedges or hash(u) <= hash(v) or \
                       v not in data or u not in data[v]:
                        self.add_edge(u,v)
        elif format == 'elliptic_curve_congruence':
            from sage.rings.arith import lcm, prime_divisors
            from sage.rings.fast_arith import prime_range
            from sage.misc.misc import prod
            for i in xrange(self.order()):
                for j in xrange(i):
                    E = curves[i]
                    F = curves[j]
                    M = E.conductor()
                    N = F.conductor()
                    MN = lcm(M, N)
                    p_MN = prime_divisors(MN)
                    lim = prod([(j^(MN.ord(j)) + j^(MN.ord(j)-1)) for j in p_MN])
                    a_E = E.anlist(lim)
                    a_F = F.anlist(lim)
                    l_list = [p for p in prime_range(lim) if p not in p_MN ]
                    p_edges = l_list
                    for l in l_list:
                        n = a_E[l] - a_F[l]
                        if n != 0:
                            P = prime_divisors(n)
                            p_edges = [p for p in p_edges if p in P]
                    if len(p_edges) > 0:
                        self.add_edge(E.cremona_label(), F.cremona_label(), str(p_edges)[1:-1])
        else:
            assert format == 'int'
        self._pos = pos
        self._boundary = boundary
        name = kwds.get('name', None)
        if format != 'Graph' or name is not None:
            self.name(name)

    ### Formats

    def graph6_string(self):
        """
        Returns the graph6 representation of the graph as an ASCII string.
        Only valid for simple (no loops, multiple edges) graphs on 0 to
        262143 vertices.

        EXAMPLES::

            sage: G = graphs.KrackhardtKiteGraph()
            sage: G.graph6_string()
            'IvUqwK@?G'
        """
        n = self.order()
        if n > 262143:
            raise ValueError('graph6 format supports graphs on 0 to 262143 vertices only.')
        elif self.has_loops() or self.has_multiple_edges():
            raise ValueError('graph6 format supports only simple graphs (no loops, no multiple edges)')
        else:
            return generic_graph_pyx.N(n) + generic_graph_pyx.R(self._bit_vector())

    def sparse6_string(self):
        """
        Returns the sparse6 representation of the graph as an ASCII string.
        Only valid for undirected graphs on 0 to 262143 vertices, but loops
        and multiple edges are permitted.

        EXAMPLES::

            sage: G = graphs.BullGraph()
            sage: G.sparse6_string()
            ':Da@en'

        ::

            sage: G = Graph()
            sage: G.sparse6_string()
            ':?'

        ::

            sage: G = Graph(loops=True, multiedges=True,sparse=True)
            sage: Graph(':?',sparse=True) == G
            True
        """
        n = self.order()
        if n == 0:
            return ':?'
        if n > 262143:
            raise ValueError('sparse6 format supports graphs on 0 to 262143 vertices only.')
        else:
            vertices = self.vertices()
            n = len(vertices)
            edges = self.edges(labels=False)
            for i in range(len(edges)): # replace edge labels with natural numbers (by index in vertices)
                edges[i] = (vertices.index(edges[i][0]),vertices.index(edges[i][1]))
            # order edges
            edges.sort(compare_edges)

            # encode bit vector
            from math import ceil
            from sage.misc.functional import log
            k = int(ceil(log(n,2)))
            v = 0
            i = 0
            m = 0
            s = ''
            while m < len(edges):
                if edges[m][1] > v + 1:
                    sp = generic_graph_pyx.binary(edges[m][1])
                    sp = '0'*(k-len(sp)) + sp
                    s += '1' + sp
                    v = edges[m][1]
                elif edges[m][1] == v + 1:
                    sp = generic_graph_pyx.binary(edges[m][0])
                    sp = '0'*(k-len(sp)) + sp
                    s += '1' + sp
                    v += 1
                    m += 1
                else:
                    sp = generic_graph_pyx.binary(edges[m][0])
                    sp = '0'*(k-len(sp)) + sp
                    s += '0' + sp
                    m += 1

            # encode s as a 6-string, as in R(x), but padding with 1's
            # pad on the right to make a multiple of 6
            s = s + ( '1' * ((6 - len(s))%6) )

            # split into groups of 6, and convert numbers to decimal, adding 63
            six_bits = ''
            for i in range(len(s)/6):
                six_bits += chr( int( s[6*i:6*(i+1)], 2) + 63 )
            return ':' + generic_graph_pyx.N(n) + six_bits

    ### Attributes

    def is_directed(self):
        """
        Since graph is undirected, returns False.

        EXAMPLES::

            sage: Graph().is_directed()
            False
        """
        return False

    ### Properties

    def eulerian_circuit(self, return_vertices=False, labels=True):
        """
        Return a list of edges forming an eulerian circuit if one exists.
        Otherwise return False.

        This is implemented using Fleury's algorithm. This could be
        extended to find eulerian paths too (check for existence and make
        sure you start on an odd-degree vertex if one exists).

        INPUT:


        -  ``return_vertices`` - optionally provide a list of
           vertices for the path

        -  ``labels`` - whether to return edges with labels
           (3-tuples)


        OUTPUT: either ([edges], [vertices]) or [edges] of an Eulerian
        circuit

        EXAMPLES::

            sage: g=graphs.CycleGraph(5);
            sage: g.eulerian_circuit()
            [(0, 1, None), (1, 2, None), (2, 3, None), (3, 4, None), (4, 0, None)]
            sage: g.eulerian_circuit(labels=False)
            [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
            sage: g = graphs.CompleteGraph(7)
            sage: edges, vertices = g.eulerian_circuit(return_vertices=True)
            sage: vertices
            [0, 1, 2, 0, 3, 1, 4, 0, 5, 1, 6, 2, 3, 4, 2, 5, 3, 6, 4, 5, 6, 0]
            sage: graphs.CompleteGraph(4).eulerian_circuit()
            False
        """
        if not self.is_eulerian():
            return False

        edge_list = []
        vertex_list = []
        from copy import copy
        g = copy(self)

        # Get first vertex
        v = g.vertex_iterator().next()
        vertex_list.append(v)
        while g.size()>0:
            for e in g.edges_incident(v, labels=labels):
                g.delete_edge(e)
                if g.is_connected():
                    break
                else:
                    g.add_edge(e)
            else:
                # Our only choice is a cut edge
                g.delete_edge(e)
                g.delete_vertex(v)
            # the following code is here so that we don't rely the
            # order of vertices in the edge tuple.
            if v == e[0]:
                v = e[1]
            else:
                v = e[0]
                e = (e[1], e[0]) + e[2:]
            edge_list.append(e)
            vertex_list.append(v)

        if return_vertices:
            return edge_list, vertex_list
        else:
            return edge_list

    def is_bipartite(self):
        """
        Returns True if graph G is bipartite, False if not.

        Traverse the graph G with depth-first-search and color nodes. This
        function uses the corresponding NetworkX function.

        EXAMPLES::

            sage: graphs.CycleGraph(4).is_bipartite()
            True
            sage: graphs.CycleGraph(5).is_bipartite()
            False
        """
        try:
            self.bipartite_color()
            return True
        except:
            return False

    def degree_constrained_subgraph(self, bounds=None):
        r"""
        Returns a degree-constrained subgraph.

        Given a graph `G` and two functions `f, g:V(G)\rightarrow \mathbb Z`
        such that `f \leq g`, a degree-constrained subgraph in `G` is
        a subgraph `G' \subseteq G` such that for any vertex `v \in G`,
        `f(v) \leq d_{G'}(v) \leq g(v)`.

        INPUT:

        - ``bounds`` -- Two possibilities :
            - A dictionary whose keys are the vertices, and values a pair of
              real values ``(min,max)`` corresponding to the values `(f(v),g(v))`.
            - A function associating to each vertex a pair of
              real values ``(min,max)`` corresponding to the values `(f(v),g(v))`.

        OUTPUT:

        - When a solution exists, this method outputs the degree-constained subgraph
          as a Graph object.
        - When no solution exists, returns ``False``

        NOTES:

        - This algorithm computes the degree-constrained subgraph of minimum weight.
        - If the graph's edges are weighted, these are taken into account.
        - This problem can be solved in polynomial time.

        EXAMPLES:

        Is there a perfect matching in an even cycle? ::

            sage: g = graphs.CycleGraph(6)
            sage: bounds = lambda x: [1,1]
            sage: m = g.degree_constrained_subgraph(bounds=bounds) # optional - requires GLPK or CBC
            sage: m.size() #optional
            3
        """

        from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException

        p = MixedIntegerLinearProgram(maximization=False)
        b = p.new_variable()

        reorder = lambda x: (min(x[0],x[1]),max(x[0],x[1]),x[2])

        if bounds is None:
            raise ValueError,"The `bounds` keyword can not be equal to None"
        elif isinstance(bounds,dict):
            f_bounds = lambda x: bounds[x]
        else:
            f_bounds = bounds


        if self.weighted():
            weight = lambda x: x[2] if x[2] is not None else 1
        else:
            weight = lambda x: 1

        for v in self:
            minimum,maximum = f_bounds(v)
            p.add_constraint(sum([ b[reorder(e)]*weight(e) for e in self.edges_incident(v)]),min=minimum, max=maximum)

        p.set_objective(sum([ b[reorder(e)]*weight(e) for e in self.edge_iterator()]))
        p.set_binary(b)

        try:
            p.solve()
            g = self.copy()
            b = p.get_values(b)
            g.delete_edges([e for e in g.edge_iterator() if b[reorder(e)] == 0])
            return g


        except MIPSolverException:
            return False





    ### Orientations

    def strong_orientation(self):
        r"""
        Returns a strongly connected orientation of the current graph.
        ( cf. http://en.wikipedia.org/wiki/Strongly_connected_component )

        An orientation of a an undirected graph is a digraph obtained by
        giving an unique direction to each of its edges. An orientation
        is said to be strong if there is a directed path between each
        pair of vertices.

        If the graph is 2-edge-connected, a strongly connected orientation
        can be found in linear time. If the given graph is not 2-connected,
        the orientation returned will ensure that each 2-connected component
        has a strongly connected orientation.

        OUTPUT:

        A digraph representing an orientation of the current graph.

        NOTES:

        - This method assumes the graph is connected.
        - This algorithm works in O(m).

        EXAMPLE:

        For a 2-regular graph, a strong orientation gives to each vertex
        an out-degree equal to 1::

            sage: g = graphs.CycleGraph(5)
            sage: g.strong_orientation().out_degree()
            [1, 1, 1, 1, 1]

        The Petersen Graph is 2-edge connected. It then has a strongly
        connected orientation::

            sage: g = graphs.PetersenGraph()
            sage: o = g.strong_orientation()
            sage: len(o.strongly_connected_components())
            1

        The same goes for the CubeGraph in any dimension ::

            sage: for i in range(2,5):
            ...     g = graphs.CubeGraph(i)
            ...     o = g.strong_orientation()
            ...     len(o.strongly_connected_components())
            1
            1
            1

        """
        from sage.graphs.all import DiGraph
        d = DiGraph(multiedges=self.allows_multiple_edges())

        id = {}
        i = 0

        # The algorithm works through a depth-first search. Any edge
        # used in the depth-first search is oriented in the direction
        # in which it has been used. All the other edges are oriented
        # backward

        v = self.vertex_iterator().next()
        seen = {}
        i=1

        # Time at which the vertices have been discovered
        seen[v] = i

        # indicates the stack of edges to explore
        next = self.edges_incident(v)

        while next:
            e = next.pop(-1)
            # We assume e[0] to be a `seen` vertex
            e = e if seen.get(e[0],False) != False else (e[1],e[0],e[2])

            # If we discovered a new vertex
            if seen.get(e[1],False) == False:
                d.add_edge(e)
                next.extend([ee for ee in self.edges_incident(e[1]) if (((e[0],e[1]) != (ee[0],ee[1])) and ((e[0],e[1]) != (ee[1],ee[0])))])
                i+=1
                seen[e[1]]=i

            # Else, we orient the edges backward
            else:
                if seen[e[0]] < seen[e[1]]:
                    d.add_edge((e[1],e[0],e[2]))
                else:
                    d.add_edge(e)

        # Case of multiple edges. If another edge has already been inserted, we add the new one
        # in the opposite direction.
        tmp = None
        for e in self.multiple_edges():
            if tmp == (e[0],e[1]):
                if d.has_edge(e[0].e[1]):
                    d.add_edge(e[1],e[0],e[2])
                else:
                    d.add_edge(e)
            tmp = (e[0],e[1])

        return d

    ### Coloring

    def bipartite_color(self):
        """
        Returns a dictionary with vertices as the keys and the color class
        as the values. Fails with an error if the graph is not bipartite.

        EXAMPLES::

            sage: graphs.CycleGraph(4).bipartite_color()
            {0: 1, 1: 0, 2: 1, 3: 0}
            sage: graphs.CycleGraph(5).bipartite_color()
            Traceback (most recent call last):
            ...
            RuntimeError: Graph is not bipartite.
        """
        # Straight from the NetworkX source:
        color = {}
        for u in self:
            if u in color:
                continue
            queue = [u]
            color[u] = 1
            while queue:
                v = queue.pop()
                c = 1-color[v]
                for w in self.neighbors(v):
                    if w in color:
                        if color[w] == color[v]:
                            raise RuntimeError("Graph is not bipartite.")
                    else:
                        color[w] = c
                        queue.append(w)
        return color

    def bipartite_sets(self):
        """
        Returns (X,Y) where X and Y are the nodes in each bipartite set of
        graph G. Fails with an error if graph is not bipartite.

        EXAMPLES::

            sage: graphs.CycleGraph(4).bipartite_sets()
            ([0, 2], [1, 3])
            sage: graphs.CycleGraph(5).bipartite_sets()
            Traceback (most recent call last):
            ...
            RuntimeError: Graph is not bipartite.
        """
        color = self.bipartite_color()
        left = [v for v in color if color[v] == 1]
        right = [v for v in color if color[v] == 0]
        return (left, right)

    def chromatic_polynomial(self):
        """
        Returns the chromatic polynomial of the graph G.

        EXAMPLES::

            sage: G = Graph({0:[1,2,3],1:[2]})
            sage: factor(G.chromatic_polynomial())
            (x - 2) * x * (x - 1)^2

        ::

            sage: g = graphs.trees(5).next()
            sage: g.chromatic_polynomial().factor()
            x * (x - 1)^4

        ::

            sage: seven_acre_wood = sum(graphs.trees(7), Graph())
            sage: seven_acre_wood.chromatic_polynomial()
            x^77 - 66*x^76 ... - 2515943049305400*x^60 ... - 66*x^12 + x^11

        ::

            sage: for i in range(2,7):
            ...     graphs.CompleteGraph(i).chromatic_polynomial().factor()
            (x - 1) * x
            (x - 2) * (x - 1) * x
            (x - 3) * (x - 2) * (x - 1) * x
            (x - 4) * (x - 3) * (x - 2) * (x - 1) * x
            (x - 5) * (x - 4) * (x - 3) * (x - 2) * (x - 1) * x
        """
        from sage.graphs.chrompoly import chromatic_polynomial
        return chromatic_polynomial(self)

    def chromatic_number(self, algorithm="DLX"):
        r"""
        Returns the minimal number of colors needed to color the vertices
        of the graph `G`.

        INPUT:

        - ``algorithm`` -- Select an algorithm from the following supported
          algorithms:

          - If ``algorithm="DLX"`` (default), the chromatic number is
            computed using the dancing link algorithm. It is
            inefficient speedwise to compute the chromatic number through
            the dancing link algorithm because this algorithm computes
            *all* the possible colorings to check that one exists.

          - If ``algorithm="CP"``, the chromatic number is computed
            using the coefficients of the chromatic polynomial. Again, this
            method is inefficient in terms of speed and it only useful for
            small graphs.

          - If ``algorithm="MILP"``, the chromatic number is computed
            using a mixed integer linear program. This method requires
            you to install an optional Sage package like GLPK or
            COIN-OR's CBC. Of the methods "DLX", "CP", and "MILP", the last
            method is the fastest method of the three.

        .. SEEALSO::

            For more functions related to graph coloring, see the
            module :mod:`sage.graphs.graph_coloring`.

        EXAMPLES::

            sage: G = Graph({0: [1, 2, 3], 1: [2]})
            sage: G.chromatic_number(algorithm="DLX")
            3
            sage: G.chromatic_number(algorithm="MILP") # optional - requires GLPK or CBC
            3
            sage: G.chromatic_number(algorithm="CP")
            3

        TESTS::

            sage: G = Graph({0: [1, 2, 3], 1: [2]})
            sage: G.chromatic_number(algorithm="foo")
            Traceback (most recent call last):
            ...
            ValueError: The 'algorithm' keyword must be set to either 'DLX', 'MILP' or 'CP'.
        """
        # default built-in algorithm; bad performance
        if algorithm == "DLX":
            from sage.graphs.graph_coloring import chromatic_number
            return chromatic_number(self)
        # Algorithm with good performance, but requires an optional
        # package: choose any of GLPK or CBC.
        elif algorithm == "MILP":
            from sage.graphs.graph_coloring import vertex_coloring
            return vertex_coloring(self, value_only=True)
        # another algorithm with bad performance; only good for small graphs
        elif algorithm == "CP":
            f = self.chromatic_polynomial()
            i = 0
            while f(i) == 0:
                i += 1
            return i
        else:
            raise ValueError("The 'algorithm' keyword must be set to either 'DLX', 'MILP' or 'CP'.")

    def coloring(self, algorithm="DLX", hex_colors=False):
        r"""
        Returns the first (optimal) proper vertex-coloring found.

        INPUT:

        - ``algorithm`` -- Select an algorithm from the following supported
          algorithms:

          - If ``algorithm="DLX"`` (default), the chromatic number is computed
            using the dancing link algorithm.

          - If ``algorithm="MILP"``, the chromatic number is computed using
            a mixed integer linear program. This algorithm requires you to
            install an optional Sage package like GLPK or COIN-OR's CBC.

        - ``hex_colors`` -- (default: ``False``) if ``True``, return a
          dictionary which can easily be used for plotting.

        .. SEEALSO::

            For more functions related to graph coloring, see the
            module :mod:`sage.graphs.graph_coloring`.

        EXAMPLES::

            sage: G = Graph("Fooba")
            sage: P = G.coloring(algorithm="MILP"); P  # optional - requires GLPK or CBC
            [[2, 1, 3], [0, 6, 5], [4]]
            sage: P = G.coloring(algorithm="DLX"); P
            [[1, 2, 3], [0, 5, 6], [4]]
            sage: G.plot(partition=P)
            sage: H = G.coloring(hex_colors=True, algorithm="MILP") # optional - requires GLPK or CBC
            sage: for c in sorted(H.keys()):                        # optional - requires GLPK or CBC
            ...       print c, H[c]                                 # optional - requires GLPK or CBC
            #0000ff [4]
            #00ff00 [0, 6, 5]
            #ff0000 [2, 1, 3]
            sage: H = G.coloring(hex_colors=True, algorithm="DLX")
            sage: for c in sorted(H.keys()):
            ...       print c, H[c]
            #0000ff [4]
            #00ff00 [1, 2, 3]
            #ff0000 [0, 5, 6]
            sage: G.plot(vertex_colors=H)

        TESTS::

            sage: G.coloring(algorithm="foo")
            Traceback (most recent call last):
            ...
            ValueError: The 'algorithm' keyword must be set to either 'DLX' or 'MILP'.
        """
        if algorithm == "MILP":
            from sage.graphs.graph_coloring import vertex_coloring
            return vertex_coloring(self, hex_colors=hex_colors)
        elif algorithm == "DLX":
            from sage.graphs.graph_coloring import first_coloring
            return first_coloring(self, hex_colors=hex_colors)
        else:
            raise ValueError("The 'algorithm' keyword must be set to either 'DLX' or 'MILP'.")

    def independent_set_of_representatives(self,family):
        r"""
        Returns an independent set of representatives.

        Given a graph `G` and and a family `F=\{F_i:i\in [1,...,k]\}` of subsets
        of ``g.vertices()``, an ISR ( Independent Set of Reprersentatives ) is an
        assignation of a vertex `v_i\in F_i` to each set `F_i` such that
        `v_i != v_j` if `i<j` ( they are represdentatives ) and the set
        `\cup_{i}v_i` is an independent set in `G`.

        It generalizes, for example, graph coloring and graph list coloring.

        ( See [AhaBerZiv07]_ for more informations )

        INPUT :

        - ``family`` -- A list of lists defining the family `F`
          ( actually, a Family of subsets of ``G.vertices()`` )

        OUTPUT :

        - A list whose `i^{\mbox{th}}` element is the representativeof the `i^{\mbox{th}}`
          element of the ``family`` list. If there is no ISR, ``None`` is returned

        EXAMPLE :

        For a bipartite graph missing one edge, the solution is as expected::

           sage: g = graphs.CompleteBipartiteGraph(3,3)
           sage: g.delete_edge(1,4)
           sage: g.independent_set_of_representatives([[0,1,2],[3,4,5]]) # optional - requires GLPK or CBC
           [1, 4]

        The Petersen Graph is 3-colorable, which can be expressed as an
        independent set of representatives problem : take 3 disjoint copies
        of the Petersen Graph, each one representing one color. Then take
        as a partition of the set of vertices the family defined by the three
        copies of each vertex. The ISR of such a family
        defines a 3-coloring ::

            sage: g = 3 * graphs.PetersenGraph()
            sage: n = g.order()/3
            sage: f = [[i,i+n,i+2*n] for i in xrange(n)]
            sage: isr = g.independent_set_of_representatives(f)   # optional - requires GLPK or CBC
            sage: c = [floor(i/n) for i in isr]                   # optional - requires GLPK or CBC
            sage: color_classes = [[],[],[]]                      # optional - requires GLPK or CBC
            sage: for v,i in enumerate(c):                        # optional - requires GLPK or CBC
            ...     color_classes[i].append(v)                    # optional - requires GLPK or CBC
            sage: for classs in color_classes:                    # optional - requires GLPK or CBC
            ...     g.subgraph(classs).size() == 0                # optional - requires GLPK or CBC
            True
            True
            True

        REFERENCE:

        .. [AhaBerZiv07] R. Aharoni and E. Berger and R. Ziv
          Independent systems of representatives in weighted graphs
          Combinatorica vol 27, num 3, p253--267
          2007

        """

        from sage.numerical.mip import MixedIntegerLinearProgram
        p=MixedIntegerLinearProgram()

        # Boolean variable indicating whether the vertex
        # is the representative of some set
        vertex_taken=p.new_variable()

        # Boolean variable in two dimension whose first
        # element is a vertex and whose second element
        # is one of the sets given as arguments.
        # When true, indicated that the vertex is the representent
        # of the corresponding set

        classss=p.new_variable(dim=2)

        # Associates to the vertices the classes
        # to which they belong

        lists=dict([(v,[]) for v in self.vertex_iterator()])
        for i,f in enumerate(family):
            [lists[v].append(i) for v in f]

            # a classss has exactly one representant
            p.add_constraint(sum([classss[v][i] for v in f]),max=1,min=1)

        # A vertex represents at most one classss (vertex_taken is binary), and
        # vertex_taken[v]==1 if v is the representative of some classss

        [p.add_constraint(sum([classss[v][i] for i in lists[v]])-vertex_taken[v],max=0) for v in self.vertex_iterator()]

        # Two adjacent vertices can not both be representants of a set

        for (u,v) in self.edges(labels=None):
            p.add_constraint(vertex_taken[u]+vertex_taken[v],max=1)

        p.set_objective(None)

        p.set_binary(vertex_taken)
        p.set_binary(classss)

        try:
            p.solve()
        except:
            return None

        classss=p.get_values(classss)

        repr=[]
        for i,f in enumerate(family):
            for v in f:
                if classss[v][i]==1:
                    repr.append(v)
                    break

        return repr

    ### Centrality

    def centrality_betweenness(self, normalized=True):
        r"""
        Returns the betweenness centrality (fraction of number of shortest
        paths that go through each vertex) as a dictionary keyed by
        vertices. The betweenness is normalized by default to be in range
        (0,1). This wraps NetworkX's implementation of the algorithm
        described in [Brandes2003]_.

        Measures of the centrality of a vertex within a graph determine the
        relative importance of that vertex to its graph. Vertices that
        occur on more shortest paths between other vertices have higher
        betweenness than vertices that occur on less.

        INPUT:


        -  ``normalized`` - boolean (default True) - if set to
           False, result is not normalized.


        REFERENCE:

        .. [Brandes2003] Ulrik Brandes. (2003). Faster Evaluation of
          Shortest-Path Based Centrality Indices. [Online] Available:
          http://citeseer.nj.nec.com/brandes00faster.html

        EXAMPLES::

            sage: (graphs.ChvatalGraph()).centrality_betweenness()
            {0: 0.069696969696969688, 1: 0.069696969696969688, 2: 0.060606060606060601, 3: 0.060606060606060601, 4: 0.069696969696969688, 5: 0.069696969696969688, 6: 0.060606060606060601, 7: 0.060606060606060601, 8: 0.060606060606060601, 9: 0.060606060606060601, 10: 0.060606060606060601, 11: 0.060606060606060601}
            sage: (graphs.ChvatalGraph()).centrality_betweenness(normalized=False)
            {0: 7.6666666666666661, 1: 7.6666666666666661, 2: 6.6666666666666661, 3: 6.6666666666666661, 4: 7.6666666666666661, 5: 7.6666666666666661, 6: 6.6666666666666661, 7: 6.6666666666666661, 8: 6.6666666666666661, 9: 6.6666666666666661, 10: 6.6666666666666661, 11: 6.6666666666666661}
            sage: D = DiGraph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: D.show(figsize=[2,2])
            sage: D = D.to_undirected()
            sage: D.show(figsize=[2,2])
            sage: D.centrality_betweenness()
            {0: 0.16666666666666666, 1: 0.16666666666666666, 2: 0.0, 3: 0.0}
        """
        import networkx
        return networkx.betweenness_centrality(self.networkx_graph(copy=False), normalized)

    def centrality_degree(self, v=None):
        r"""
        Returns the degree centrality (fraction of vertices connected to)
        as a dictionary of values keyed by vertex. The degree centrality is
        normalized to be in range (0,1).

        Measures of the centrality of a vertex within a graph determine the
        relative importance of that vertex to its graph. Degree centrality
        measures the number of links incident upon a vertex.

        INPUT:


        -  ``v`` - a vertex label (to find degree centrality of
           only one vertex)


        EXAMPLES::

            sage: (graphs.ChvatalGraph()).centrality_degree()
            {0: 0.36363636363636365, 1: 0.36363636363636365, 2: 0.36363636363636365, 3: 0.36363636363636365, 4: 0.36363636363636365, 5: 0.36363636363636365, 6: 0.36363636363636365, 7: 0.36363636363636365, 8: 0.36363636363636365, 9: 0.36363636363636365, 10: 0.36363636363636365, 11: 0.36363636363636365}
            sage: D = DiGraph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: D.show(figsize=[2,2])
            sage: D = D.to_undirected()
            sage: D.show(figsize=[2,2])
            sage: D.centrality_degree()
            {0: 1.0, 1: 1.0, 2: 0.66666666666666663, 3: 0.66666666666666663}
            sage: D.centrality_degree(v=1)
            1.0
        """
        import networkx
        return networkx.degree_centrality(self.networkx_graph(copy=False), v)

    def centrality_closeness(self, v=None):
        r"""
        Returns the closeness centrality (1/average distance to all
        vertices) as a dictionary of values keyed by vertex. The degree
        centrality is normalized to be in range (0,1).

        Measures of the centrality of a vertex within a graph determine the
        relative importance of that vertex to its graph. 'Closeness
        centrality may be defined as the total graph-theoretic distance of
        a given vertex from all other vertices... Closeness is an inverse
        measure of centrality in that a larger value indicates a less
        central actor while a smaller value indicates a more central
        actor,' [Borgatti95]_.

        INPUT:


        -  ``v`` - a vertex label (to find degree centrality of
           only one vertex)


        REFERENCE:

        .. [Borgatti95] Stephen P Borgatti. (1995). Centrality and AIDS.
          [Online] Available:
          http://www.analytictech.com/networks/centaids.htm

        EXAMPLES::

            sage: (graphs.ChvatalGraph()).centrality_closeness()
            {0: 0.61111111111111116, 1: 0.61111111111111116, 2: 0.61111111111111116, 3: 0.61111111111111116, 4: 0.61111111111111116, 5: 0.61111111111111116, 6: 0.61111111111111116, 7: 0.61111111111111116, 8: 0.61111111111111116, 9: 0.61111111111111116, 10: 0.61111111111111116, 11: 0.61111111111111116}
            sage: D = DiGraph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: D.show(figsize=[2,2])
            sage: D = D.to_undirected()
            sage: D.show(figsize=[2,2])
            sage: D.centrality_closeness()
            {0: 1.0, 1: 1.0, 2: 0.75, 3: 0.75}
            sage: D.centrality_closeness(v=1)
            1.0
        """
        import networkx
        return networkx.closeness_centrality(self.networkx_graph(copy=False), v)

    ### Constructors

    def to_directed(self, implementation='c_graph', sparse=None):
        """
        Returns a directed version of the graph. A single edge becomes two
        edges, one in each direction.

        EXAMPLES::

            sage: graphs.PetersenGraph().to_directed()
            Petersen graph: Digraph on 10 vertices
        """
        if sparse is None:
            from sage.graphs.base.dense_graph import DenseGraphBackend
            sparse = (not isinstance(self._backend, DenseGraphBackend))
        from sage.graphs.all import DiGraph
        D = DiGraph(name=self.name(), pos=self._pos, boundary=self._boundary,
                    multiedges=self.allows_multiple_edges(),
                    implementation=implementation, sparse=sparse)
        D.name(self.name())
        D.add_vertices(self.vertex_iterator())
        for u,v,l in self.edge_iterator():
            D.add_edge(u,v,l)
            D.add_edge(v,u,l)
        if hasattr(self, '_embedding'):
            from copy import copy
            D._embedding = copy(self._embedding)
        D._weighted = self._weighted
        return D

    def to_undirected(self):
        """
        Since the graph is already undirected, simply returns a copy of
        itself.

        EXAMPLES::

            sage: graphs.PetersenGraph().to_undirected()
            Petersen graph: Graph on 10 vertices
        """
        from copy import copy
        return copy(self)

    ### Visualization

    def write_to_eps(self, filename, iterations=50):
        r"""
        Writes a plot of the graph to filename in eps format.

        It is relatively simple to include this file in a latex document:

        INPUT: filename


        -  ``iterations`` - how many iterations of the spring
           layout algorithm to go through, if applicable


        ``usepackagegraphics`` must appear before the beginning of the
        document, and ``includegraphics filename.eps`` will include it in your latex
        doc. Note: you cannot use pdflatex to print the resulting document,
        use TeX and Ghostscript or something similar instead.

        EXAMPLES::

            sage: P = graphs.PetersenGraph()
            sage: P.write_to_eps(tmp_dir() + 'sage.eps')
        """
        from sage.graphs.print_graphs import print_graph_eps
        if self._pos is None:
            pos = generic_graph_pyx.spring_layout_fast(self, iterations=iterations)
        else:
            pos = self._pos
            keys = pos.keys()
            for v in self.vertices():
                if v not in keys:
                    pos = generic_graph_pyx.spring_layout_fast(self, iterations=iterations)
                    break
        xmin = 0.0
        ymin = 0.0
        xmax = -1.0
        ymax = -1.0
        for v in pos:
            x,y = pos[v]
            if (x > xmax):
                xmax = x
            if (x < xmin):
                xmin = x
            if (y > ymax):
                ymax = y
            if (y < ymin):
                ymin = y
        for v in pos:
            pos[v][0] = 1.8*(pos[v][0] - xmin)/(xmax - xmin) - 0.9
            pos[v][1] = 1.8*(pos[v][1] - ymin)/(ymax - ymin) - 0.9
        if filename[-4:] != '.eps':
            filename += '.eps'
        f = open(filename, 'w')
        f.write( print_graph_eps(self.vertices(), self.edge_iterator(), pos) )
        f.close()

    def graphviz_string(self):
       r"""
       Returns a representation in the DOT language, ready to render in
       graphviz.

       REFERENCES:

       - http://www.graphviz.org/doc/info/lang.html

       EXAMPLES::

           sage: G = Graph({0:{1:None,2:None}, 1:{0:None,2:None}, 2:{0:None,1:None,3:'foo'}, 3:{2:'foo'}},sparse=True)
           sage: s = G.graphviz_string()
           sage: s
           'graph {\n"0";"1";"2";"3";\n"0"--"1";"0"--"2";"1"--"2";"2"--"3"[label="foo"];\n}'
       """
       return self._graphviz_string_helper("graph", "--") # edge_string is "--" for undirected graphs

    ### Cliques

    def cliques_maximal(self):
        """
        Returns the list of all maximal cliques, with each clique represented
        by a list of vertices. A clique is an induced complete subgraph, and a
        maximal clique is one not contained in a larger one.

        NOTES:

         - Currently only implemented for undirected graphs. Use to_undirected
           to convert a digraph to an undirected graph.

        ALGORITHM:

        This function is based on NetworkX's implementation of the Bron and
        Kerbosch Algorithm [BroKer1973]_.

        REFERENCE:

        .. [BroKer1973] Coen Bron and Joep Kerbosch. (1973). Algorithm 457:
          Finding All Cliques of an Undirected Graph. Commun. ACM. v
          16. n 9.  pages 575-577. ACM Press. [Online] Available:
          http://www.ram.org/computing/rambin/rambin.html

        EXAMPLES::

            sage: graphs.ChvatalGraph().cliques_maximal()
            [[0, 1], [0, 4], [0, 6], [0, 9], [2, 1], [2, 3], [2, 6], [2, 8], [3, 4], [3, 7], [3, 9], [5, 1], [5, 4], [5, 10], [5, 11], [7, 1], [7, 8], [7, 11], [8, 4], [8, 10], [10, 6], [10, 9], [11, 6], [11, 9]]
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_maximal()
            [[0, 1, 2], [0, 1, 3]]
            sage: C=graphs.PetersenGraph()
            sage: C.cliques_maximal()
            [[0, 1], [0, 4], [0, 5], [2, 1], [2, 3], [2, 7], [3, 4], [3, 8], [6, 1], [6, 8], [6, 9], [7, 5], [7, 9], [8, 5], [9, 4]]
            sage: C = Graph('DJ{')
            sage: C.cliques_maximal()
            [[4, 1, 2, 3], [4, 0]]

        """
        import networkx.cliques
        return networkx.cliques.find_cliques(self.networkx_graph(copy=False))

    def cliques(self):
        """
        (Deprecated) alias for ``cliques_maximal``. See that function for more
        details.

        EXAMPLE::

            sage: C = Graph('DJ{')
            sage: C.cliques()
            doctest:...: DeprecationWarning: The function 'cliques' has been deprecated. Use 'cliques_maximal' or 'cliques_maximum'.
            [[4, 1, 2, 3], [4, 0]]

        """
        from sage.misc.misc import deprecation
        deprecation("The function 'cliques' has been deprecated. Use " + \
                    "'cliques_maximal' or 'cliques_maximum'.")
        return self.cliques_maximal()

    def cliques_maximum(self):
        """
        Returns the list of all maximum cliques, with each clique represented
        by a list of vertices. A clique is an induced complete subgraph, and a
        maximum clique is one of maximal order.

        NOTES:

        - Currently only implemented for undirected graphs. Use to_undirected
          to convert a digraph to an undirected graph.

        ALGORITHM:

        This function is based on Cliquer [NisOst2003]_.

        REFERENCE:

        .. [NisOst2003] Sampo Niskanen and Patric R. J. Ostergard,
          "Cliquer User's Guide, Version 1.0," Communications Laboratory,
          Helsinki University of Technology, Espoo, Finland,
          Tech. Rep. T48, 2003.

        EXAMPLES::

            sage: graphs.ChvatalGraph().cliques_maximum()
            [[0, 1], [0, 4], [0, 6], [0, 9], [1, 2], [1, 5], [1, 7], [2, 3], [2, 6], [2, 8], [3, 4], [3, 7], [3, 9], [4, 5], [4, 8], [5, 10], [5, 11], [6, 10], [6, 11], [7, 8], [7, 11], [8, 10], [9, 10], [9, 11]]
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_maximum()
            [[0, 1, 2], [0, 1, 3]]
            sage: C=graphs.PetersenGraph()
            sage: C.cliques_maximum()
            [[0, 1], [0, 4], [0, 5], [1, 2], [1, 6], [2, 3], [2, 7], [3, 4], [3, 8], [4, 9], [5, 7], [5, 8], [6, 8], [6, 9], [7, 9]]
            sage: C = Graph('DJ{')
            sage: C.cliques_maximum()
            [[1, 2, 3, 4]]

        """
        from sage.graphs.cliquer import all_max_clique
        return sorted(all_max_clique(self))

    def clique_maximum(self):
        """
        Returns the vertex set of a maximal order complete subgraph.

        NOTE:

         - Currently only implemented for undirected graphs. Use to_undirected
           to convert a digraph to an undirected graph.

        ALGORITHM:

        This function is based on Cliquer [NisOst2003]_.

        EXAMPLES::

            sage: C=graphs.PetersenGraph()
            sage: C.clique_maximum()
            [7, 9]
            sage: C = Graph('DJ{')
            sage: C.clique_maximum()
            [1, 2, 3, 4]

        """
        from sage.graphs.cliquer import max_clique
        return max_clique(self)

    def clique_number(self, algorithm="cliquer", cliques=None):
        r"""
        Returns the order of the largest clique of the graph (the clique
        number).

        NOTE:

         - Currently only implemented for undirected graphs. Use ``to_undirected``
           to convert a digraph to an undirected graph.

        INPUT:

         - ``algorithm`` - either ``cliquer`` or ``networkx``

           - ``cliquer`` - This wraps the C program Cliquer [NisOst2003]_.

           - ``networkx`` - This function is based on NetworkX's implementation
             of the Bron and Kerbosch Algorithm [BroKer1973]_.

         - ``cliques`` - an optional list of cliques that can be input if
           already computed. Ignored unless ``algorithm=='networkx'``.

        ALGORITHM:

        This function is based on Cliquer [NisOst2003]_ and [BroKer1973]_.

        EXAMPLES::

            sage: C = Graph('DJ{')
            sage: C.clique_number()
            4
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.clique_number()
            3
        """
        if algorithm=="cliquer":
            from sage.graphs.cliquer import clique_number
            return clique_number(self)
        elif algorithm=="networkx":
            import networkx.cliques
            return networkx.cliques.graph_clique_number(self.networkx_graph(copy=False),cliques)
        else:
            raise NotImplementedError("Only 'networkx' and 'cliquer' are supported.")

    def cliques_number_of(self, vertices=None, cliques=None, with_labels=False):
        """
        Returns a list of the number of maximal cliques containing each
        vertex. (Returns a single value if only one input vertex).

        NOTES:

         - Currently only implemented for undirected graphs. Use to_undirected
           to convert a digraph to an undirected graph.

        INPUT:

        -  ``vertices`` - the vertices to inspect (default is
           entire graph)

        -  ``with_labels`` - (boolean) default False returns
           list as above True returns a dictionary keyed by vertex labels

        -  ``cliques`` - list of cliques (if already
           computed)


        EXAMPLES::

            sage: C = Graph('DJ{')
            sage: C.cliques_number_of()
            [1, 1, 1, 1, 2]
            sage: E = C.cliques_maximal()
            sage: E
            [[4, 1, 2, 3], [4, 0]]
            sage: C.cliques_number_of(cliques=E)
            [1, 1, 1, 1, 2]
            sage: F = graphs.Grid2dGraph(2,3)
            sage: X = F.cliques_number_of(with_labels=True)
            sage: for v in sorted(X.iterkeys()):
            ...    print v, X[v]
            (0, 0) 2
            (0, 1) 3
            (0, 2) 2
            (1, 0) 2
            (1, 1) 3
            (1, 2) 2
            sage: F.cliques_number_of(vertices=[(0, 1), (1, 2)])
            [3, 2]
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_number_of()
            [2, 2, 1, 1]
        """
        import networkx.cliques
        return networkx.cliques.number_of_cliques(self.networkx_graph(copy=False), vertices, cliques, with_labels)

    def cliques_get_max_clique_graph(self, name=''):
        """
        Returns a graph constructed with maximal cliques as vertices, and
        edges between maximal cliques with common members in the original
        graph.

        NOTES:

         - Currently only implemented for undirected graphs. Use to_undirected
           to convert a digraph to an undirected graph.

        INPUT:

        -  ``name`` - The name of the new graph.

        EXAMPLES::

            sage: (graphs.ChvatalGraph()).cliques_get_max_clique_graph()
            Graph on 24 vertices
            sage: ((graphs.ChvatalGraph()).cliques_get_max_clique_graph()).show(figsize=[2,2], vertex_size=20, vertex_labels=False)
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_get_max_clique_graph()
            Graph on 2 vertices
            sage: (G.cliques_get_max_clique_graph()).show(figsize=[2,2])
        """
        import networkx.cliques
        return Graph(networkx.cliques.make_max_clique_graph(self.networkx_graph(copy=False), name=name, create_using=networkx.xgraph.XGraph()))

    def cliques_get_clique_bipartite(self, **kwds):
        """
        Returns a bipartite graph constructed such that maximal cliques are the
        right vertices and the left vertices are retained from the given
        graph. Right and left vertices are connected if the bottom vertex
        belongs to the clique represented by a top vertex.

        NOTES:

         - Currently only implemented for undirected graphs. Use to_undirected
           to convert a digraph to an undirected graph.

        EXAMPLES::

            sage: (graphs.ChvatalGraph()).cliques_get_clique_bipartite()
            Bipartite graph on 36 vertices
            sage: ((graphs.ChvatalGraph()).cliques_get_clique_bipartite()).show(figsize=[2,2], vertex_size=20, vertex_labels=False)
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_get_clique_bipartite()
            Bipartite graph on 6 vertices
            sage: (G.cliques_get_clique_bipartite()).show(figsize=[2,2])
        """
        import networkx.cliques
        from bipartite_graph import BipartiteGraph
        return BipartiteGraph(networkx.cliques.make_clique_bipartite(self.networkx_graph(copy=False), **kwds))

    def independent_set(self):
        """
        Returns a maximal independent set, which is a set of vertices which
        induces an empty subgraph. Uses Cliquer [NisOst2003]_.

        NOTES:

         - Currently only implemented for undirected graphs. Use to_undirected
           to convert a digraph to an undirected graph.

        EXAMPLES::

            sage: C=graphs.PetersenGraph()
            sage: C.independent_set()
            [0, 3, 6, 7]
        """
        from sage.graphs.cliquer import max_clique
        return max_clique(self.complement())

    def cliques_vertex_clique_number(self, algorithm="cliquer", vertices=None,
                                     with_labels=False, cliques=None):
        r"""
        Returns a list of sizes of the largest maximal cliques containing
        each vertex. (Returns a single value if only one input vertex).

        NOTES:

         - Currently only implemented for undirected graphs. Use to_undirected
           to convert a digraph to an undirected graph.

        INPUT:

         - ``algorithm`` - either ``cliquer`` or ``networkx``

           - ``cliquer`` - This wraps the C program Cliquer [NisOst2003]_.

           - ``networkx`` - This function is based on NetworkX's implementation
                of the Bron and Kerbosch Algorithm [BroKer1973]_.

        -  ``vertices`` - the vertices to inspect (default is entire graph).
           Ignored unless ``algorithm=='networkx'``.

        -  ``with_labels`` - (boolean) default False returns list as above
           True returns a dictionary keyed by vertex labels. Ignored unless
           ``algorithm=='networkx'``.

        -  ``cliques`` - list of cliques (if already computed).  Ignored unless
           ``algorithm=='networkx'``.

        EXAMPLES::

            sage: C = Graph('DJ{')
            sage: C.cliques_vertex_clique_number()
            [2, 4, 4, 4, 4]
            sage: E = C.cliques_maximal()
            sage: E
            [[4, 1, 2, 3], [4, 0]]
            sage: C.cliques_vertex_clique_number(cliques=E,algorithm="networkx")
            [2, 4, 4, 4, 4]
            sage: F = graphs.Grid2dGraph(2,3)
            sage: X = F.cliques_vertex_clique_number(with_labels=True,algorithm="networkx")
            sage: for v in sorted(X.iterkeys()):
            ...    print v, X[v]
            (0, 0) 2
            (0, 1) 2
            (0, 2) 2
            (1, 0) 2
            (1, 1) 2
            (1, 2) 2
            sage: F.cliques_vertex_clique_number(vertices=[(0, 1), (1, 2)])
            [2, 2]
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_vertex_clique_number()
            [3, 3, 3, 3]

        """

        if algorithm=="cliquer":
            from sage.graphs.cliquer import clique_number
            if vertices==None:
                vertices=self
            value=[]
            for v in vertices:
                value.append(1+clique_number(self.subgraph(self.neighbors(v))))
                self.subgraph(self.neighbors(v)).plot()
            return value
        elif algorithm=="networkx":
            import networkx.cliques
            return networkx.cliques.node_clique_number(self.networkx_graph(copy=False), vertices, with_labels, cliques)
        else:
            raise NotImplementedError("Only 'networkx' and 'cliquer' are supported.")

    def cliques_containing_vertex(self, vertices=None, cliques=None, with_labels=False):
        """
        Returns the cliques containing each vertex, represented as a list
        of lists. (Returns a single list if only one input vertex).

        NOTE:

         - Currently only implemented for undirected graphs. Use to_undirected
           to convert a digraph to an undirected graph.

        INPUT:

        -  ``vertices`` - the vertices to inspect (default is
           entire graph)

        -  ``with_labels`` - (boolean) default False returns
           list as above True returns a dictionary keyed by vertex labels

        -  ``cliques`` - list of cliques (if already
           computed)

        EXAMPLES::

            sage: C = Graph('DJ{')
            sage: C.cliques_containing_vertex()
            [[[4, 0]], [[4, 1, 2, 3]], [[4, 1, 2, 3]], [[4, 1, 2, 3]], [[4, 1, 2, 3], [4, 0]]]
            sage: E = C.cliques_maximal()
            sage: E
            [[4, 1, 2, 3], [4, 0]]
            sage: C.cliques_containing_vertex(cliques=E)
            [[[4, 0]], [[4, 1, 2, 3]], [[4, 1, 2, 3]], [[4, 1, 2, 3]], [[4, 1, 2, 3], [4, 0]]]
            sage: F = graphs.Grid2dGraph(2,3)
            sage: X = F.cliques_containing_vertex(with_labels=True)
            sage: for v in sorted(X.iterkeys()):
            ...    print v, X[v]
            (0, 0) [[(0, 1), (0, 0)], [(1, 0), (0, 0)]]
            (0, 1) [[(0, 1), (0, 0)], [(0, 1), (1, 1)], [(0, 1), (0, 2)]]
            (0, 2) [[(0, 1), (0, 2)], [(1, 2), (0, 2)]]
            (1, 0) [[(1, 0), (0, 0)], [(1, 0), (1, 1)]]
            (1, 1) [[(0, 1), (1, 1)], [(1, 2), (1, 1)], [(1, 0), (1, 1)]]
            (1, 2) [[(1, 2), (1, 1)], [(1, 2), (0, 2)]]
            sage: F.cliques_containing_vertex(vertices=[(0, 1), (1, 2)])
            [[[(0, 1), (0, 0)], [(0, 1), (1, 1)], [(0, 1), (0, 2)]], [[(1, 2), (1, 1)], [(1, 2), (0, 2)]]]
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_containing_vertex()
            [[[0, 1, 2], [0, 1, 3]], [[0, 1, 2], [0, 1, 3]], [[0, 1, 2]], [[0, 1, 3]]]

        """
        import networkx.cliques
        return networkx.cliques.cliques_containing_node(self.networkx_graph(copy=False), vertices, cliques, with_labels)

    def clique_complex(self):
        """
        Returns the clique complex of self. This is the largest simplicial complex on
        the vertices of self whose 1-skeleton is self.

        This is only makes sense for undirected simple graphs.

        EXAMPLES::

            sage: g = Graph({0:[1,2],1:[2],4:[]})
            sage: g.clique_complex()
            Simplicial complex with vertex set (0, 1, 2, 4) and facets {(4,), (0, 1, 2)}

            sage: h = Graph({0:[1,2,3,4],1:[2,3,4],2:[3]})
            sage: x = h.clique_complex()
            sage: x
            Simplicial complex with vertex set (0, 1, 2, 3, 4) and facets {(0, 1, 4), (0, 1, 2, 3)}
            sage: i = x.graph()
            sage: i==h
            True
            sage: x==i.clique_complex()
            True

        """
        if self.is_directed() or self.has_loops() or self.has_multiple_edges():
            raise ValueError("Self must be an undirected simple graph to have a clique complex.")
        import sage.homology.simplicial_complex
        C = sage.homology.simplicial_complex.SimplicialComplex(self.vertices(),self.cliques_maximal(),maximality_check=True)
        C._graph = self
        return C

    ### Miscellaneous

    def min_spanning_tree(self, weight_function=lambda e: 1,
                          algorithm='Kruskal',
                          starting_vertex=None ):
        """
        Returns the edges of a minimum spanning tree, if one exists,
        otherwise returns False.

        INPUT:


        -  ``weight_function`` - A function that takes an edge
           and returns a numeric weight. Defaults to assigning each edge a
           weight of 1.

        -  ``algorithm`` - Three variants of algorithms are
           implemented: 'Kruskal', 'Prim fringe', and 'Prim edge' (the last
           two are variants of Prim's algorithm). Defaults to 'Kruskal'.
           Currently, 'Prim fringe' ignores the labels on the edges.

        -  ``starting_vertex`` - The vertex with which to
           start Prim's algorithm.


        OUTPUT: the edges of a minimum spanning tree.

        EXAMPLES::

            sage: g=graphs.CompleteGraph(5)
            sage: len(g.min_spanning_tree())
            4
            sage: weight = lambda e: 1/( (e[0]+1)*(e[1]+1) )
            sage: g.min_spanning_tree(weight_function=weight)
            [(3, 4, None), (2, 4, None), (1, 4, None), (0, 4, None)]
            sage: g.min_spanning_tree(algorithm='Prim edge', starting_vertex=2, weight_function=weight)
            [(2, 4, None), (3, 4, None), (1, 3, None), (0, 4, None)]
            sage: g.min_spanning_tree(algorithm='Prim fringe', starting_vertex=2, weight_function=weight)
            [(4, 2), (3, 4), (1, 4), (0, 4)]
        """
        if self.is_connected()==False:
            return False

        if algorithm=='Kruskal':
            # Kruskal's algorithm
            edges=[]
            sorted_edges_iterator=iter(sorted(self.edges(), key=weight_function))
            union_find = dict([(v,None) for v in self.vertex_iterator()])
            while len(edges) < self.order()-1:
                # get next edge
                e=sorted_edges_iterator.next()
                components=[]
                for start_v in e[0:2]:
                    v=start_v
                    children=[]

                    # Find the component a vertex lives in.
                    while union_find[v] != None:
                        children.append(v)
                        v=union_find[v]

                    # Compress the paths as much as we can for
                    # efficiency reasons.
                    for child in children:
                        union_find[child]=v

                    components.append(v)

                if components[0]!=components[1]:
                    # put in edge
                    edges.append(e)

                    # Union the components by making one the parent of the
                    # other.
                    union_find[components[0]]=components[1]

            return edges

        elif algorithm=='Prim fringe':
            if starting_vertex is None:
                v = self.vertex_iterator().next()
            else:
                v = starting_vertex
            tree=set([v])
            edges=[]

            # initialize fringe_list with v's neighbors.  fringe_list
            # contains fringe_vertex: (vertex_in_tree, weight) for each
            # fringe vertex
            fringe_list=dict([u,(v,weight_function((v,u)))] for u in self[v])

            for i in xrange(self.order()-1):
                # Find the smallest-weight fringe vertex
                v=min(fringe_list,key=lambda x: fringe_list[x][1])
                edges.append((v,fringe_list[v][0]))
                tree.add(v)
                fringe_list.pop(v)

                # Update fringe list
                for neighbor in [u for u in self[v] if u not in tree]:
                    w=weight_function((v,neighbor))
                    if neighbor not in fringe_list or \
                           (neighbor in fringe_list and fringe_list[neighbor][1]>w):
                        fringe_list[neighbor]=(v,weight_function((v,neighbor)))
            return edges

        elif algorithm=='Prim edge':
            if starting_vertex is None:
                v = self.vertex_iterator().next()
            else:
                v = starting_vertex
            sorted_edges=sorted(self.edges(), key=weight_function)
            tree=set([v])
            edges=[]

            for i in xrange(self.order()-1):
                # Find a minimum-weight edge connecting a vertex in
                # the tree to something outside the tree.  Remove the
                # edges between tree vertices for efficiency.

                for i in xrange(len(sorted_edges)):
                    e=sorted_edges[i]
                    v0,v1=e[0],e[1]
                    if v0 in tree:
                        if v1 not in tree:
                            edges.append(e)
                            sorted_edges[i:i+1]=[]
                            tree.add(v1)
                            break
                        else:
                            sorted_edges[i:i+1]=[]
                    elif v1 in tree:
                        edges.append(e)
                        sorted_edges[i:i+1]=[]
                        tree.add(v0)
                        break
            return edges
        else:
            raise NotImplementedError, "Minimum Spanning Tree algorithm '%s' is not implemented."%algorithm

    def two_factor_petersen(self):
        r"""
        Returns a decomposition of the graph into 2-factors.

        Petersen's 2-factor decomposition theorem asserts that any
        `2r`-regular graph `G` can be decomposed into 2-factors.
        Equivalently, it means that the edges of any `2r`-regular
        graphs can be partitionned in `r` sets `C_1,\dots,C_r` such
        that for all `i`, the set `C_i` is a disjoint union of cycles
        ( a 2-regular graph ).

        As any graph of maximal degree `\Delta` can be completed into
        a regular graph of degree `2\lceil\frac\Delta 2\rceil`, this
        result also means that the edges of any graph of degree `\Delta`
        can be partitionned in `r=2\lceil\frac\Delta 2\rceil` sets
        `C_1,\dots,C_r` such that for all `i`, the set `C_i` is a
        graph of maximal degree `2` ( a disjoint union of paths
        and cycles ).

        EXAMPLE:

        The Complete Graph on `7` vertices is a `6`-regular graph, so it can
        be edge-partitionned into `2`-regular graphs::

            sage: g = graphs.CompleteGraph(7)
            sage: classes = g.two_factor_petersen()  # optional - requires GLPK or CBC
            sage: for c in classes:                  # optional - requires GLPK or CBC
            ...     gg = Graph()                     # optional - requires GLPK or CBC
            ...     gg.add_edges(c)                  # optional - requires GLPK or CBC
            ...     print max(gg.degree())<=2        # optional - requires GLPK or CBC
            True
            True
            True
            sage: Set(set(classes[0]) | set(classes[1]) | set(classes[2])).cardinality() == g.size() # optional - requires GLPK or CBC
            True

        ::

            sage: g = graphs.CirculantGraph(24, [7, 11])
            sage: cl = g.two_factor_petersen()                     # optional - GLPK or CBC
            sage: g.plot(edge_colors={'black':cl[0], 'red':cl[1]}) # optional - GLPK or CBC

        """

        d = self.eulerian_orientation()

        # This new graph is bipartite, and built the following way :
        #
        # To each vertex v of the digraph are associated two vertices,
        # a sink (-1,v) and a source (1,v)
        # Any edge (u,v) in the digraph is then added as ((-1,u),(1,v))

        from sage.graphs.graph import Graph
        g = Graph()
        g.add_edges([((-1,u),(1,v)) for (u,v) in d.edge_iterator(labels=None)])

        # This new bipartite graph is now edge_colored
        from sage.graphs.graph_coloring import edge_coloring
        classes = edge_coloring(g,log=1)

        # The edges in the classes are of the form ((-1,u),(1,v))
        # and have to be translated back to (u,v)
        classes_b = []
        for c in classes:
            classes_b.append([(u,v) for ((uu,u),(vv,v),t) in c])

        return classes_b


def compare_edges(x, y):
    """
    Compare edge x to edge y, return -1 if x y, 1 if x y, else 0.

    EXAMPLES::

        sage: G = graphs.PetersenGraph()
        sage: E = G.edges()
        sage: from sage.graphs.graph import compare_edges
        sage: compare_edges(E[0], E[2])
        -1
        sage: compare_edges(E[0], E[1])
        -1
        sage: compare_edges(E[0], E[0])
        0
        sage: compare_edges(E[1], E[0])
        1
    """
    if x[1] < y[1]:
        return -1
    elif x[1] > y[1]:
        return 1
    elif x[1] == y[1]:
        if x[0] < y[0]:
            return -1
        if x[0] > y[0]:
            return 1
        else:
            return 0

