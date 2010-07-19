r"""
Undirected graphs

This module implements functions and operations involving undirected
graphs.

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

- Nicolas M. Thiery (2010-02): graph layout code refactoring, dot2tex/graphviz interface


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

The NetworkX graph is essentially a dictionary of dictionaries of dictionaries::

    sage: N.adj
    {0: {1: {}, 4: {}, 5: {}}, 1: {0: {}, 2: {}, 6: {}}, 2: {1: {}, 3: {}, 7: {}}, 3: {8: {}, 2: {}, 4: {}}, 4: {0: {}, 9: {}, 3: {}}, 5: {0: {}, 8: {}, 7: {}}, 6: {8: {}, 1: {}, 9: {}}, 7: {9: {}, 2: {}, 5: {}}, 8: {3: {}, 5: {}, 6: {}}, 9: {4: {}, 6: {}, 7: {}}}

Each dictionary key is a vertex label, and each key in the
following dictionary is a neighbor of that vertex. In undirected
graphs, there is redundancy: for example, the dictionary containing
the entry ``1: {2: {}}`` implies it must contain
``{2: {1: {}}``. The innermost entry of ``{}`` is
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

-  a list of edges::

       sage: g = Graph([(1,3),(3,8),(5,2)])
       sage: g
       Graph on 5 vertices

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

    A graph is a set of vertices connected by edges. See also the
    `Wikipedia article on graphs <http://en.wikipedia.org/wiki/Graph_(mathematics)>`_.

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

    #. a list of edges, or labelled edges::

          sage: g = Graph([(1,3),(3,8),(5,2)])
          sage: g
          Graph on 5 vertices

          ::

          sage: g = Graph([(1,2,"Peace"),(7,-9,"and"),(77,2, "Love")])
          sage: g
          Graph on 5 vertices

    #. A NetworkX MultiGraph::

          sage: import networkx
          sage: g = networkx.MultiGraph({0:[1,2,3], 2:[4]})
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

        The positions are copied when the graph is built from
        another graph ::

            sage: g = graphs.PetersenGraph()
            sage: h = Graph(g)
            sage: g.get_pos() == h.get_pos()
            True

        Or from a DiGraph ::

            sage: d = DiGraph(g)
            sage: h = Graph(d)
            sage: g.get_pos() == h.get_pos()
            True

        Invalid sequence of edges given as an input (they do not all
        have the same length)::

            sage: g = Graph([(1,2),(2,3),(2,3,4)])
            Traceback (most recent call last):
            ...
            ValueError: Edges input must all follow the same format.


        Two different labels given for the same edge in a graph
        without multiple edges::

            sage: g = Graph([(1,2,3),(1,2,4)], multiedges = False)
            Traceback (most recent call last):
            ...
            ValueError: Two different labels given for the same edge in a graph without multiple edges.
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
            if isinstance(data, (networkx.DiGraph, networkx.MultiDiGraph)):
                data = data.to_undirected()
                format = 'NX'
            elif isinstance(data, (networkx.Graph, networkx.MultiGraph)):
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

                # The edges are not labelled
                if len(data[0]) == 2:
                    data = {}
                    for u,v in edges:
                        if not u in data:
                            data[u] = []
                        if not v in data:
                            data[v] = []
                        data[u].append(v)
                        data[v].append(u)

                    format = 'dict_of_lists'

                # The edges are labelled
                elif len(data[0]) == 3:
                    data = {}
                    for u,v,l in edges:

                        if not u in data:
                            data[u] = {}
                        if not v in data:
                            data[v] = {}

                        # Now the keys exists, and are dictionaries ...


                        # If we notice for the first time that there
                        # are multiple edges, we update the whole
                        # dictionary so that data[u][v] is a list

                        if (multiedges is None and
                            (u in data[v])):
                            multiedges = True
                            for uu, dd in data.iteritems():
                                for vv, ddd in dd.iteritems():
                                    dd[vv] = [ddd]

                        # If multiedges is set to False while the same
                        # edge is present in the list with different
                        # values of its label
                        elif (multiedges is False and
                            (u in data[v] and data[v][u] != l)):
                                raise ValueError("MULTIEDGE")

                        # Now we are behaving as if multiedges == None
                        # means multiedges = False. If something bad
                        # happens later, the whole dictionary will be
                        # updated anyway

                        if multiedges is True:
                            if v not in data[u]:
                                data[u][v] = []
                                data[v][u] = []

                            data[u][v].append(l)
                            data[v][u].append(l)

                        else:
                            data[u][v] = l
                            data[v][u] = l

                    format = 'dict_of_dicts'

            except ValueError as ve:
                if str(ve) == "MULTIEDGE":
                    raise ValueError("Two different labels given for the same edge in a graph without multiple edges.")
                else:
                    raise ValueError("Edges input must all follow the same format.")


        if format is None:
            import networkx
            data = networkx.MultiGraph(data)
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
                        multiedges = True
                    if loops is None:
                        loops = True
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
                if multiedges:
                    self._backend = NetworkXGraphBackend(networkx.MultiGraph())
                else:
                    self._backend = NetworkXGraphBackend(networkx.Graph())
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

    def is_even_hole_free(self, certificate = False):
        r"""
        Tests whether ``self`` contains an induced even hole.

        A Hole is a cycle of length at least 4 (included). It is said
        to be even (resp. odd) if its length is even (resp. odd).

        Even-hole-free graphs always contain a bisimplicial vertex,
        which ensures that their chromatic number is at most twice
        their clique number [ABCHRS08]_.

        INPUT:

        - ``certificate`` (boolean) -- When ``certificate = False``,
          this method only returns ``True`` or ``False``. If
          ``certificate = True``, the subgraph found is returned
          instead of ``False``.

        EXAMPLE:

        Is the Petersen Graph even-hole-free ::

            sage: g = graphs.PetersenGraph()
            sage: g.is_even_hole_free()
            False

        As any chordal graph is hole-free, interval graphs behave the
        same way::

            sage: g = graphs.RandomInterval(20)
            sage: g.is_even_hole_free()
            True

        It is clear, though, that a random Bipartite Graph which is
        not a forest has an even hole::

            sage: g = graphs.RandomBipartite(10, 10, .5)
            sage: g.is_even_hole_free() and not g.is_forest()
            False

        We can check the certificate returned is indeed an even
        cycle::

            sage: if not g.is_forest():
            ...      cycle = g.is_even_hole_free(certificate = True)
            ...      if cycle.order() % 2 == 1:
            ...          print "Error !"
            ...      if not cycle.is_isomorphic(
            ...             graphs.CycleGraph(cycle.order())):
            ...          print "Error !"
            ...
            sage: print "Everything is Fine !"
            Everything is Fine !

        TESTS:

        Bug reported in #9925, and fixed by #9420::

            sage: g = Graph(':SiBFGaCEF_@CE`DEGH`CEFGaCDGaCDEHaDEF`CEH`ABCDEF')
            sage: g.is_even_hole_free()
            False
            sage: g.is_even_hole_free(certificate = True)
            Subgraph of (): Looped multi-graph on 4 vertices

        Making sure there are no other counter-examples around ::

            sage: t = lambda x : (Graph(x).is_forest() or
            ...         isinstance(Graph(x).is_even_hole_free(certificate = True),Graph))
            sage: all( t(graphs.RandomBipartite(10,10,.5)) for i in range(100) )
            True

        REFERENCE:

        .. [ABCHRS08] L. Addario-Berry, M. Chudnovsky, F. Havet, B. Reed, P. Seymour
          Bisimplicial  vertices in even-hole-free graphs
          Journal of Combinatorial Theory, Series B
          vol 98, n.6 pp 1119-1164, 2008
        """
        from sage.graphs.graph_generators import GraphGenerators

        girth = self.girth()

        if girth > self.order():
            start = 4

        elif girth % 2 == 0:
            if not certificate:
                return False
            start = girth

        else:
            start = girth + 1

        while start <= self.order():


            subgraph = self.subgraph_search(GraphGenerators().CycleGraph(start), induced = True)

            if not subgraph is None:
                if certificate:
                    return subgraph
                else:
                    return False

            start = start + 2

        return True

    def is_odd_hole_free(self, certificate = False):
        r"""
        Tests whether ``self`` contains an induced odd hole.

        A Hole is a cycle of length at least 4 (included). It is said
        to be even (resp. odd) if its length is even (resp. odd).

        It is interesting to notice that while it is polynomial to
        check whether a graph has an odd hole or an odd antihole [CRST06]_, it is
        not known whether testing for one of these two cases
        independently is polynomial too.

        INPUT:

        - ``certificate`` (boolean) -- When ``certificate = False``,
          this method only returns ``True`` or ``False``. If
          ``certificate = True``, the subgraph found is returned
          instead of ``False``.

        EXAMPLE:

        Is the Petersen Graph odd-hole-free ::

            sage: g = graphs.PetersenGraph()
            sage: g.is_odd_hole_free()
            False

        Which was to be expected, as its girth is 5 ::

            sage: g.girth()
            5

        We can check the certificate returned is indeed a 5-cycle::

            sage: cycle = g.is_odd_hole_free(certificate = True)
            sage: cycle.is_isomorphic(graphs.CycleGraph(5))
            True

        As any chordal graph is hole-free, no interval graph has an odd hole::

            sage: g = graphs.RandomInterval(20)
            sage: g.is_odd_hole_free()
            True

        REFERENCES:

        .. [CRST06] M. Chudnovsky, G. Cornuejols, X. Liu, P. Seymour, K. Vuskovic
          Recognizing berge graphs
          Combinatorica vol 25, n 2, pages 143--186
          2005
        """
        from sage.graphs.graph_generators import GraphGenerators

        girth = self.girth()

        if girth > self.order() or girth == 3:
            start = 5

        elif girth % 2 == 1:
            if not certificate:
                return False
            start = girth

        else:
            start = girth + 1

        while start <= self.order():

            subgraph = self.subgraph_search(GraphGenerators().CycleGraph(start), induced = True)

            if not subgraph is None:
                if certificate:
                    return subgraph
                else:
                    return False

            start = start + 2

        return True

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
            [(0, 1, {}), (1, 2, {}), (2, 3, {}), (3, 4, {}), (4, 0, {})]
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

    def is_line_graph(self, certificate = False):
        r"""
        Tests wether the graph is a line graph.

        INPUT:

        - ``certificate`` (boolean) -- whether to return a certificate
          when the graph is *not* a line graph. When ``certificate``
          is set to ``True``, and if the graph is not a line graph,
          the method returns a subgraph isomorphic to one of the 9
          forbidden induced subgraphs of a line graph (instead of the
          usual ``False``)


        TODO:

        This methods sequentially tests each of the forbidden
        subgraphs, which is a very slow method. There exist much
        better algorithms, including those which are actually able to
        return a graph whose line graph is isomorphic to the given
        graph.

        EXAMPLES:

        A complete graph is always the line graph of a star::

            sage: graphs.CompleteGraph(5).is_line_graph()
            True

        The Petersen Graph not being claw-free, it is not a line
        graph:

            sage: graphs.PetersenGraph().is_line_graph()
            False

        This is indeed the subgraph returned::

            sage: C = graphs.PetersenGraph().is_line_graph(certificate = True)
            sage: C.is_isomorphic(graphs.ClawGraph())
            True
        """
        from sage.graphs.graph_generators import graphs

        for g in graphs.line_graph_forbidden_subgraphs():
            h = self.subgraph_search(g, induced = True)
            if h is not None:
                if certificate:
                    return h
                else:
                    return False

        return True

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

    def is_triangle_free(self):
        r"""
        Returns whether ``self`` is triangle-free

        EXAMPLE:

        The Petersen Graph is triangle-free::

            sage: g = graphs.PetersenGraph()
            sage: g.is_triangle_free()
            True

        or a complete Bipartite Graph::

            sage: g = graphs.CompleteBipartiteGraph(5,6)
            sage: g.is_triangle_free()
            True

        a tripartite graph, though, contains many triangles ::

            sage: g = (3 * graphs.CompleteGraph(5)).complement()
            sage: g.is_triangle_free()
            False
        """

        from sage.graphs.graph_generators import graphs

        return (self.subgraph_search(graphs.CompleteGraph(3)) is None)

    def is_split(self):
        r"""
        Returns ``True`` if the graph is a Split graph, ``False`` otherwise.

        A Graph `G` is said to be a split graph if its vertices `V(G)`
        can be partitioned into two sets `K` and `I` such that the
        vertices of `K` induce a complete graphe, and those of `I` are
        an independent set.

        There is a simple test to check whether a graph is a split
        graph (see, for instance, the book "Graph Classes, a survey"
        [GraphClasses]_ page 203) :

        Given the degree sequence `d_1 \geq ... \geq d_n` of `G`, a graph
        is a split graph if and only if :

        .. MATH::

            \sum_{i=1}^\omega d_i = \omega (\omega - 1) + \sum_{i=\omega + 1}^nd_i

        where `\omega = max \{i:d_i\geq i-1\}`.


        EXAMPLES:

        Split graphs are, in particular, chordal graphs. Hence, The Petersen graph
        can not be split::

            sage: graphs.PetersenGraph().is_split()
            False

        We can easily build some "random" split graph by creating a
        complete graph, and adding vertices only connected
        to some random vertices of the clique::

            sage: g = graphs.CompleteGraph(10)
            sage: sets = Subsets(Set(range(10)))
            sage: for i in range(10, 25):
            ...      g.add_edges([(i,k) for k in sets.random_element()])
            sage: g.is_split()
            True

        REFERENCES:

        .. [GraphClasses] A. Brandstadt, VB Le and JP Spinrad
          Graph classes: a survey
          SIAM Monographs on Discrete Mathematics and Applications},
          1999
        """

        # our degree sequence is numbered from 0 to n-1, so to avoid
        # any mistake, let's fix it :-)
        degree_sequence = [0] + sorted(self.degree(), reverse = True)

        for (i, d) in enumerate(degree_sequence):
            if d >= i - 1:
                omega = i
            else:
                break

        left = sum(degree_sequence[:omega + 1])
        right = omega * (omega - 1) + sum(degree_sequence[omega + 1:])

        return left == right


    def is_perfect(self, certificate = False):
        r"""
        Tests whether the graph is perfect.

        A graph `G` is said to be perfect if `\chi(H)=\omega(H)` hold
        for any induced subgraph `H\subseteq_i G` (and so for `G`
        itself, too), where `\chi(H)` represents the chromatic number
        of `H`, and `\omega(H)` its clique number. The Strong Perfect
        Graph Theorem [SPGT]_ gives another characterization of
        perfect graphs:

        A graph is perfect if and only if it contains no odd hole
        (cycle on an odd number `k` of vertices, `k>3`) nor any odd
        antihole (complement of a hole) as an induced subgraph.

        INPUT:

        - ``certificate`` (boolean) -- whether to return
          a certificate (default : ``False``)

        OUTPUT:

        When ``certificate = False``, this function returns
        a boolean value. When ``certificate = True``, it returns
        a subgraph of ``self`` isomorphic to an odd hole or an odd
        antihole if any, and ``None`` otherwise.

        EXAMPLE:

        A Bipartite Graph is always perfect ::

            sage: g = graphs.RandomBipartite(8,4,.5)
            sage: g.is_perfect()
            True

        Interval Graphs, which are chordal graphs, too ::

            sage: g =  graphs.RandomInterval(7)
            sage: g.is_perfect()
            True

        The PetersenGraph, which is triangle-free and
        has chromatic number 3 is obviously not perfect::

            sage: g = graphs.PetersenGraph()
            sage: g.is_perfect()
            False

        We can obtain an induced 5-cycle as a certificate::

            sage: g.is_perfect(certificate = True)
            Subgraph of (Petersen graph): Graph on 5 vertices

        REFERENCES:

        .. [SPGT] M. Chudnovsky, N. Robertson, P. Seymour, R. Thomas.
          The strong perfect graph theorem
          Annals of Mathematics
          vol 164, number 1, pages 51--230
          2006
        """

        from sage.graphs.bipartite_graph import BipartiteGraph

        if isinstance(self, BipartiteGraph) or self.is_bipartite():
            return True if not certificate else None

        self_complement = self.complement()

        start_complement = self_complement.girth()
        start_self = self.girth()

        from sage.graphs.graph_generators import graphs


        # In these cases, we know the graph is no perfect.
        if start_self == 5:
            if certificate:
                return self.subgraph_search(graphs.CycleGraph(5), induced = True)
            else:
                return False

        if start_complement == 5:
            if certificate:
                return self_complement.subgraph_search(graphs.CycleGraph(5), induced = True).complement()
            else:
                return False

        # We are only looking for odd holes of size at least 5
        from sage.rings.finite_rings.integer_mod import Mod

        start = lambda x : (x+1) if Mod(x,2) == 0 else ( 5 if x == 3 else x )

        # these values are possibly the infinity !!!!

        start_self = start(start_self) if start_self <= self.order() else self.order()+2
        start_complement = start(start_complement) if start_complement <= self.order() else self.order()+2

        counter_example = None

        for i in range(min(start_self, start_complement), self.order()+1,2):

            # trying in self
            if i >= start_self:
                counter_example = self.subgraph_search(graphs.CycleGraph(i), induced = True)

                if counter_example is not None:
                    break

            # trying in the complement
            if i >= start_complement:
                counter_example = self.subgraph_search(graphs.CycleGraph(i), induced = True)
                if counter_example is not None:
                    counter_example = counter_example.complement()
                    break

        if certificate:
            return counter_example
        else:
            return counter_example is None


    def degree_constrained_subgraph(self, bounds=None, solver=None, verbose=0):
        r"""
        Returns a degree-constrained subgraph.

        Given a graph `G` and two functions `f, g:V(G)\rightarrow \mathbb Z`
        such that `f \leq g`, a degree-constrained subgraph in `G` is
        a subgraph `G' \subseteq G` such that for any vertex `v \in G`,
        `f(v) \leq d_{G'}(v) \leq g(v)`.

        INPUT:

        - ``bounds`` -- (default: ``None``) Two possibilities:

          - A dictionary whose keys are the vertices, and values a pair of
            real values ``(min,max)`` corresponding to the values
            `(f(v),g(v))`.

          - A function associating to each vertex a pair of
            real values ``(min,max)`` corresponding to the values
            `(f(v),g(v))`.


        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        OUTPUT:

        - When a solution exists, this method outputs the degree-constained
          subgraph as a Graph object.

        - When no solution exists, returns ``False``.

        .. NOTE::

            - This algorithm computes the degree-constrained subgraph of minimum weight.
            - If the graph's edges are weighted, these are taken into account.
            - This problem can be solved in polynomial time.

        EXAMPLES:

        Is there a perfect matching in an even cycle? ::

            sage: g = graphs.CycleGraph(6)
            sage: bounds = lambda x: [1,1]
            sage: m = g.degree_constrained_subgraph(bounds=bounds)
            sage: m.size()
            3
        """

        from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException, Sum

        p = MixedIntegerLinearProgram(maximization=False, solver=solver)
        b = p.new_variable()

        reorder = lambda x,y: (x,y) if x<y else (y,x)

        if bounds is None:
            raise ValueError,"The `bounds` keyword can not be equal to None"
        elif isinstance(bounds,dict):
            f_bounds = lambda x: bounds[x]
        else:
            f_bounds = bounds


        if self.weighted():
            from sage.rings.real_mpfr import RR
            weight = lambda x: x if x in RR else 1
        else:
            weight = lambda x: 1

        for v in self:
            minimum,maximum = f_bounds(v)
            p.add_constraint(Sum([ b[reorder(x,y)]*weight(l) for x,y,l in self.edges_incident(v)]), min=minimum, max=maximum)

        p.set_objective(Sum([ b[reorder(x,y)]*weight(l) for x,y,l in self.edge_iterator()]))
        p.set_binary(b)

        try:
            p.solve(log=verbose)
            g = self.copy()
            b = p.get_values(b)
            g.delete_edges([(x,y) for x,y,_ in g.edge_iterator() if b[reorder(x,y)] < 0.5])
            return g


        except MIPSolverException:
            return False


    ### Orientations

    def strong_orientation(self):
        r"""
        Returns a strongly connected orientation of the current graph. See
        also the
        `Wikipedia article on strongly connected component <http://en.wikipedia.org/wiki/Strongly_connected_component>`_.

        An orientation of an undirected graph is a digraph obtained by
        giving an unique direction to each of its edges. An orientation
        is said to be strong if there is a directed path between each
        pair of vertices.

        If the graph is 2-edge-connected, a strongly connected orientation
        can be found in linear time. If the given graph is not 2-connected,
        the orientation returned will ensure that each 2-connected component
        has a strongly connected orientation.

        OUTPUT:

        A digraph representing an orientation of the current graph.

        .. NOTE::

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

    def bounded_outdegree_orientation(self, bound):
        r"""
        Computes an orientation of ``self`` such that every vertex `v`
        has out-degree less than `b(v)`

        INPUT:

        - ``bound`` -- Maximum bound on the out-degree. Can be of
          three different types :

             * An integer `k`. In this case, computes an orientation
               whose maximum out-degree is less than `k`.

             * A dictionary associating to each vertex its associated
               maximum out-degree.

             * A function associating to each vertex its associated
               maximum out-degree.

        OUTPUT:

        A DiGraph representing the orientation if it exists. A
        ``ValueError`` exception is raised otherwise.

        ALGORITHM:

        The problem is solved through a maximum flow :

        Given a graph `G`, we create a ``DiGraph`` `D` defined on
        `E(G)\cup V(G)\cup \{s,t\}`. We then link `s` to all of `V(G)`
        (these edges having a capacity equal to the bound associated
        to each element of `V(G)`), and all the elements of `E(G)` to
        `t` . We then link each `v \in V(G)` to each of its incident
        edges in `G`. A maximum integer flow of value `|E(G)|`
        corresponds to an admissible orientation of `G`. Otherwise,
        none exists.

        EXAMPLES:

        There is always an orientation of a graph `G` such that a
        vertex `v` has out-degree at most `\lceil \frac {d(v)} 2
        \rceil`::

            sage: g = graphs.RandomGNP(40, .4)
            sage: b = lambda v : ceil(g.degree(v)/2)
            sage: D = g.bounded_outdegree_orientation(b)
            sage: all( D.out_degree(v) <= b(v) for v in g )
            True


        Chvatal's graph, being 4-regular, can be oriented in such a
        way that its maximum out-degree is 2::

            sage: g = graphs.ChvatalGraph()
            sage: D = g.bounded_outdegree_orientation(2)
            sage: max(D.out_degree())
            2

        For any graph `G`, it is possible to compute an orientation
        such that the maximum out-degree is at most the maximum
        average degree of `G` divided by 2. Anything less, though, is
        impossible.

            sage: g = graphs.RandomGNP(40, .4)
            sage: mad = g.maximum_average_degree()

        Hence this is possible ::

            sage: d = g.bounded_outdegree_orientation(ceil(mad/2))

        While this is not::

            sage: try:
            ...      g.bounded_outdegree_orientation(ceil(mad/2-1))
            ...      print "Error"
            ... except ValueError:
            ...       pass

        TESTS:

        As previously for random graphs, but more intensively::

            sage: for i in xrange(30):                                   # long
            ...       g = graphs.RandomGNP(40, .4)                       # long
            ...       b = lambda v : ceil(g.degree(v)/2)                 # long
            ...       D = g.bounded_outdegree_orientation(b)             # long
            ...       if not (                                           # long
            ...            all( D.out_degree(v) <= b(v) for v in g ) or  # long
            ...            D.size() != g.size()):                        # long
            ...           print "Something wrong happened"               # long

        """
        from sage.graphs.all import DiGraph
        n = self.order()

        if n == 0:
            return DiGraph()

        vertices = self.vertices()
        vertices_id = dict(map(lambda (x,y):(y,x), list(enumerate(vertices))))

        b = {}


        # Checking the input type. We make a dictionay out of it
        if isinstance(bound, dict):
            b = bound
        else:
            try:
                b = dict(zip(vertices,map(bound, vertices)))

            except TypeError:
                b = dict(zip(vertices, [bound]*n))

        d = DiGraph()

        # Adding the edges (s,v) and ((u,v),t)
        d.add_edges( ('s', vertices_id[v], b[v]) for v in vertices)

        d.add_edges( ((vertices_id[u], vertices_id[v]), 't', 1)
                     for (u,v) in self.edges(labels=None) )

        # each v is linked to its incident edges

        for u,v in self.edges(labels = None):
            u,v = vertices_id[u], vertices_id[v]
            d.add_edge(u, (u,v), 1)
            d.add_edge(v, (u,v), 1)

        # Solving the maximum flow
        value, flow = d.flow('s','t', value_only = False, integer = True, use_edge_labels = True)

        if value != self.size():
            raise ValueError("No orientation exists for the given bound")

        D = DiGraph()
        D.add_vertices(vertices)

        # The flow graph may not contain all the vertices, if they are
        # not part of the flow...

        for u in [x for x in range(n) if x in flow]:

            for (uu,vv) in flow.neighbors_out(u):
                v = vv if vv != u else uu
                D.add_edge(vertices[u], vertices[v])

        # I do not like when a method destroys the embedding ;-)

        D.set_pos(self.get_pos())

        return D


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
            (set([0, 2]), set([1, 3]))
            sage: graphs.CycleGraph(5).bipartite_sets()
            Traceback (most recent call last):
            ...
            RuntimeError: Graph is not bipartite.
        """
        color = self.bipartite_color()
        left = set([v for v in color if color[v] == 1])
        right = set([v for v in color if color[v] == 0])
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
            sage: G.chromatic_number(algorithm="MILP")
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
            sage: P = G.coloring(algorithm="MILP"); P
            [[2, 1, 3], [0, 6, 5], [4]]
            sage: P = G.coloring(algorithm="DLX"); P
            [[1, 2, 3], [0, 5, 6], [4]]
            sage: G.plot(partition=P)
            sage: H = G.coloring(hex_colors=True, algorithm="MILP")
            sage: for c in sorted(H.keys()):
            ...       print c, H[c]
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

    def independent_set_of_representatives(self, family, solver=None, verbose=0):
        r"""
        Returns an independent set of representatives.

        Given a graph `G` and and a family `F=\{F_i:i\in [1,...,k]\}` of
        subsets of ``g.vertices()``, an Independent Set of Reprersentatives
        (ISR) is an assignation of a vertex `v_i\in F_i` to each set `F_i`
        such that `v_i != v_j` if `i<j` (they are represdentatives) and the
        set `\cup_{i}v_i` is an independent set in `G`.

        It generalizes, for example, graph coloring and graph list coloring.

        (See [AhaBerZiv07]_ for more information.)

        INPUT:

        - ``family`` -- A list of lists defining the family `F`
          (actually, a Family of subsets of ``G.vertices()``).

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        OUTPUT:

        - A list whose `i^{\mbox{th}}` element is the representativeof the
          `i^{\mbox{th}}` element of the ``family`` list. If there is no ISR,
          ``None`` is returned.

        EXAMPLES:

        For a bipartite graph missing one edge, the solution is as expected::

           sage: g = graphs.CompleteBipartiteGraph(3,3)
           sage: g.delete_edge(1,4)
           sage: g.independent_set_of_representatives([[0,1,2],[3,4,5]])
           [1, 4]

        The Petersen Graph is 3-colorable, which can be expressed as an
        independent set of representatives problem : take 3 disjoint copies
        of the Petersen Graph, each one representing one color. Then take
        as a partition of the set of vertices the family defined by the three
        copies of each vertex. The ISR of such a family
        defines a 3-coloring::

            sage: g = 3 * graphs.PetersenGraph()
            sage: n = g.order()/3
            sage: f = [[i,i+n,i+2*n] for i in xrange(n)]
            sage: isr = g.independent_set_of_representatives(f)
            sage: c = [floor(i/n) for i in isr]
            sage: color_classes = [[],[],[]]
            sage: for v,i in enumerate(c):
            ...     color_classes[i].append(v)
            sage: for classs in color_classes:
            ...     g.subgraph(classs).size() == 0
            True
            True
            True

        REFERENCE:

        .. [AhaBerZiv07] R. Aharoni and E. Berger and R. Ziv
          Independent systems of representatives in weighted graphs
          Combinatorica vol 27, num 3, p253--267
          2007

        """

        from sage.numerical.mip import MixedIntegerLinearProgram, Sum
        p=MixedIntegerLinearProgram(solver=solver)

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
            p.add_constraint(Sum([classss[v][i] for v in f]),max=1,min=1)

        # A vertex represents at most one classss (vertex_taken is binary), and
        # vertex_taken[v]==1 if v is the representative of some classss

        [p.add_constraint(Sum([classss[v][i] for i in lists[v]])-vertex_taken[v],max=0) for v in self.vertex_iterator()]

        # Two adjacent vertices can not both be representants of a set

        for (u,v) in self.edges(labels=None):
            p.add_constraint(vertex_taken[u]+vertex_taken[v],max=1)

        p.set_objective(None)

        p.set_binary(vertex_taken)
        p.set_binary(classss)

        try:
            p.solve(log=verbose)
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

    def minor(self, H, solver=None, verbose=0):
        r"""
        Returns the vertices of a minor isomorphic to `H` in the current graph.

        We say that a graph `G` has a `H`-minor (or that it has
        a graph isomorphic to `H` as a minor), if for all `h\in H`,
        there exist disjoint sets `S_h \subseteq V(G)` such that
        once the vertices of each `S_h` have been merged to create
        a new graph `G'`, this new graph contains `H` as a subgraph.

        For more information, see the
        `Wikipedia article on graph minor <http://en.wikipedia.org/wiki/Minor_%28graph_theory%29>`_.

        INPUT:

        - ``H`` -- The minor to find for in the current graph.

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        OUTPUT:

        A dictionary associating to each vertex of `H` the set of vertices
        in the current graph representing it.

        ALGORITHM:

        Mixed Integer Linear Programming

        COMPLEXITY:

        Theoretically, when `H` is fixed, testing for the existence of
        a `H`-minor is polynomial. The known algorithms are highly
        exponential in `H`, though.

        .. NOTE::

            This function can be expected to be *very* slow, especially
            where the minor does not exist.

        EXAMPLES:

        Trying to find a minor isomorphic to `K_4` in
        the `4\times 4` grid::

            sage: g = graphs.GridGraph([4,4])
            sage: h = graphs.CompleteGraph(4)
            sage: L = g.minor(h)
            sage: gg = g.subgraph(flatten(L.values(), max_level = 1))
            sage: _ = [gg.merge_vertices(l) for l in L.values() if len(l)>1]
            sage: gg.is_isomorphic(h)
            True

        We can also try to prove this way that the Petersen graph
        is not planar, as it has a `K_5` minor::

            sage: g = graphs.PetersenGraph()
            sage: K5_minor = g.minor(graphs.CompleteGraph(5))                    # long time

        And even a `K_{3,3}` minor::

            sage: K33_minor = g.minor(graphs.CompleteBipartiteGraph(3,3))        # long time

        (It is much faster to use the linear-time test of
        planarity in this situation, though.)

        As there is no cycle in a tree, looking for a `K_3` minor is useless.
        This function will raise an exception in this case::

            sage: g = graphs.RandomGNP(20,.5)
            sage: g = g.subgraph(edges = g.min_spanning_tree())
            sage: g.is_tree()
            True
            sage: L = g.minor(graphs.CompleteGraph(3))
            Traceback (most recent call last):
            ...
            ValueError: This graph has no minor isomorphic to H !
        """

        from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException, Sum
        p = MixedIntegerLinearProgram(solver=solver)

        # sorts an edge
        S = lambda (x,y) : (x,y) if x<y else (y,x)

        # rs = Representative set of a vertex
        # for h in H, v in G is such that rs[h][v] == 1 if and only if v
        # is a representant of h in self
        rs = p.new_variable(dim=2)

        for v in self:
            p.add_constraint(Sum([rs[h][v] for h in H]), max = 1)

        # We ensure that the set of representatives of a
        # vertex h contains a tree, and thus is connected

        # edges represents the edges of the tree
        edges = p.new_variable(dim = 2)

        # there can be a edge for h between two vertices
        # only if those vertices represent h
        for u,v in self.edges(labels=None):
            for h in H:
                p.add_constraint(edges[h][S((u,v))] - rs[h][u], max = 0 )
                p.add_constraint(edges[h][S((u,v))] - rs[h][v], max = 0 )

        # The number of edges of the tree in h is exactly the cardinal
        # of its representative set minus 1

        for h in H:
            p.add_constraint(Sum([edges[h][S(e)] for e in self.edges(labels=None)])-Sum([rs[h][v] for v in self]), min=-1, max=-1)

        # a tree  has no cycle
        epsilon = 1/(5*Integer(self.order()))
        r_edges = p.new_variable(dim=2)

        for h in H:
            for u,v in self.edges(labels=None):
                p.add_constraint(r_edges[h][(u,v)] + r_edges[h][(v,u)] - edges[h][S((u,v))], min = 0)

            for v in self:
                p.add_constraint(Sum([r_edges[h][(u,v)] for u in self.neighbors(v)]), max = 1-epsilon)

        # Once the representative sets are described, we must ensure
        # there are arcs corresponding to those of H between them
        h_edges = p.new_variable(dim=2)

        for h1, h2 in H.edges(labels=None):

            for v1, v2 in self.edges(labels=None):

                p.add_constraint(h_edges[(h1,h2)][S((v1,v2))] - rs[h2][v2], max = 0)
                p.add_constraint(h_edges[(h1,h2)][S((v1,v2))] - rs[h1][v1], max = 0)

                p.add_constraint(h_edges[(h2,h1)][S((v1,v2))] - rs[h1][v2], max = 0)
                p.add_constraint(h_edges[(h2,h1)][S((v1,v2))] - rs[h2][v1], max = 0)

            p.add_constraint(Sum([h_edges[(h1,h2)][S(e)] + h_edges[(h2,h1)][S(e)] for e in self.edges(labels=None) ]), min = 1)

        p.set_binary(rs)
        p.set_binary(edges)

        p.set_objective(None)

        try:
            p.solve(log=verbose)
        except MIPSolverException:
            raise ValueError("This graph has no minor isomorphic to H !")

        rs = p.get_values(rs)

        rs_dict = {}
        for h in H:
            rs_dict[h] = [v for v in self if rs[h][v]==1]

        return rs_dict

    ### Matching

    def matching_polynomial(self, complement=True, name=None):
        """
        Computes the matching polynomial of the graph `G`.

        If `p(G, k)` denotes the number of `k`-matchings (matchings with `k`
        edges) in `G`, then the matching polynomial is defined as [Godsil93]_:

        .. MATH::

            \mu(x)=\sum_{k \geq 0} (-1)^k p(G,k) x^{n-2k}


        INPUT:

        - ``complement`` - (default: ``True``) whether to use Godsil's duality
          theorem to compute the matching polynomial from that of the graphs
          complement (see ALGORITHM).

        - ``name`` - optional string for the variable name in the polynomial

        .. NOTE::

            The ``complement`` option uses matching polynomials of complete
            graphs, which are cached. So if you are crazy enough to try
            computing the matching polynomial on a graph with millions of
            vertices, you might not want to use this option, since it will end
            up caching millions of polynomials of degree in the millions.

        ALGORITHM:

        The algorithm used is a recursive one, based on the following
        observation [Godsil93]_:

        - If `e` is an edge of `G`, `G'` is the result of deleting the edge `e`,
          and `G''` is the result of deleting each vertex in `e`, then the
          matching polynomial of `G` is equal to that of `G'` minus that of
          `G''`.

          (the algorithm actually computes the *signless* matching polynomial,
          for which the recursion is the same when one replaces the substraction
          by an addition. It is then converted into the matching polynomial and
          returned)

        Depending on the value of ``complement``, Godsil's duality theorem
        [Godsil93]_ can also be used to compute `\mu(x)` :

        .. MATH::

            \mu(\overline{G}, x) = \sum_{k \geq 0} p(G,k) \mu( K_{n-2k}, x)


        Where `\overline{G}` is the complement of `G`, and `K_n` the complete
        graph on `n` vertices.

        EXAMPLES::

            sage: g = graphs.PetersenGraph()
            sage: g.matching_polynomial()
            x^10 - 15*x^8 + 75*x^6 - 145*x^4 + 90*x^2 - 6
            sage: g.matching_polynomial(complement=False)
            x^10 - 15*x^8 + 75*x^6 - 145*x^4 + 90*x^2 - 6
            sage: g.matching_polynomial(name='tom')
            tom^10 - 15*tom^8 + 75*tom^6 - 145*tom^4 + 90*tom^2 - 6
            sage: g = Graph()
            sage: L = [graphs.RandomGNP(8, .3) for i in [1..5]]
            sage: prod([h.matching_polynomial() for h in L]) == sum(L, g).matching_polynomial()
            True
        """
        from matchpoly import matching_polynomial
        return matching_polynomial(self, complement=complement, name=name)

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
            {0: 3.833333333333333, 1: 3.833333333333333, 2: 3.333333333333333, 3: 3.333333333333333, 4: 3.833333333333333, 5: 3.833333333333333, 6: 3.333333333333333, 7: 3.333333333333333, 8: 3.333333333333333, 9: 3.333333333333333, 10: 3.333333333333333, 11: 3.333333333333333}
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
        if v==None:
            return networkx.degree_centrality(self.networkx_graph(copy=False))
        else:
            return networkx.degree_centrality(self.networkx_graph(copy=False))[v]

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
            {0: 0.61111111111111..., 1: 0.61111111111111..., 2: 0.61111111111111..., 3: 0.61111111111111..., 4: 0.61111111111111..., 5: 0.61111111111111..., 6: 0.61111111111111..., 7: 0.61111111111111..., 8: 0.61111111111111..., 9: 0.61111111111111..., 10: 0.61111111111111..., 11: 0.61111111111111...}
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

    def write_to_eps(self, filename, **options):
        r"""
        Writes a plot of the graph to ``filename`` in ``eps`` format.

        INPUT:

         - ``filename`` -- a string
         - ``**options`` -- same layout options as :meth:`.layout`

        EXAMPLES::

            sage: P = graphs.PetersenGraph()
            sage: P.write_to_eps(tmp_dir() + 'sage.eps')

        It is relatively simple to include this file in a LaTeX
        document.  ``\usepackagegraphics`` must appear in the
        preamble, and ``\includegraphics{filename.eps}`` will include
        the file. To compile the document to ``pdf`` with ``pdflatex``
        the file needs first to be converted to ``pdf``, for example
        with ``ps2pdf filename.eps filename.pdf``.
        """
        from sage.graphs.print_graphs import print_graph_eps
        pos = self.layout(**options)
        [xmin, xmax, ymin, ymax] = self._layout_bounding_box(pos)
        for v in pos:
            pos[v] = (1.8*(pos[v][0] - xmin)/(xmax - xmin) - 0.9, 1.8*(pos[v][1] - ymin)/(ymax - ymin) - 0.9)
        if filename[-4:] != '.eps':
            filename += '.eps'
        f = open(filename, 'w')
        f.write( print_graph_eps(self.vertices(), self.edge_iterator(), pos) )
        f.close()


    ### Cliques

    def cliques_maximal(self):
        """
        Returns the list of all maximal cliques, with each clique represented
        by a list of vertices. A clique is an induced complete subgraph, and a
        maximal clique is one not contained in a larger one.

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
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
            [[4, 0], [4, 1, 2, 3]]

        """
        import networkx
        return sorted(networkx.find_cliques(self.networkx_graph(copy=False)))

    def cliques(self):
        """
        (Deprecated) alias for ``cliques_maximal``. See that function for more
        details.

        EXAMPLE::

            sage: C = Graph('DJ{')
            sage: C.cliques()
            doctest:...: DeprecationWarning: The function 'cliques' has been deprecated. Use 'cliques_maximal' or 'cliques_maximum'.
            [[4, 0], [4, 1, 2, 3]]

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

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
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

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
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

        .. NOTE::

            Currently only implemented for undirected graphs. Use ``to_undirected``
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
            import networkx
            return networkx.graph_clique_number(self.networkx_graph(copy=False),cliques)
        else:
            raise NotImplementedError("Only 'networkx' and 'cliquer' are supported.")

    def cliques_number_of(self, vertices=None, cliques=None):
        """
        Returns a dictionary of the number of maximal cliques containing each
        vertex, keyed by vertex. (Returns a single value if
        only one input vertex).

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
            to convert a digraph to an undirected graph.

        INPUT:

        -  ``vertices`` - the vertices to inspect (default is
           entire graph)

        -  ``cliques`` - list of cliques (if already
           computed)


        EXAMPLES::

            sage: C = Graph('DJ{')
            sage: C.cliques_number_of()
            {0: 1, 1: 1, 2: 1, 3: 1, 4: 2}
            sage: E = C.cliques_maximal()
            sage: E
            [[4, 0], [4, 1, 2, 3]]
            sage: C.cliques_number_of(cliques=E)
            {0: 1, 1: 1, 2: 1, 3: 1, 4: 2}
            sage: F = graphs.Grid2dGraph(2,3)
            sage: X = F.cliques_number_of()
            sage: for v in sorted(X.iterkeys()):
            ...    print v, X[v]
            (0, 0) 2
            (0, 1) 3
            (0, 2) 2
            (1, 0) 2
            (1, 1) 3
            (1, 2) 2
            sage: F.cliques_number_of(vertices=[(0, 1), (1, 2)])
            {(0, 1): 3, (1, 2): 2}
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_number_of()
            {0: 2, 1: 2, 2: 1, 3: 1}
        """
        import networkx
        return networkx.number_of_cliques(self.networkx_graph(copy=False), vertices, cliques)

    def cliques_get_max_clique_graph(self, name=''):
        """
        Returns a graph constructed with maximal cliques as vertices, and
        edges between maximal cliques with common members in the original
        graph.

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
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
        import networkx
        return Graph(networkx.make_max_clique_graph(self.networkx_graph(copy=False), name=name, create_using=networkx.MultiGraph()))

    def cliques_get_clique_bipartite(self, **kwds):
        """
        Returns a bipartite graph constructed such that maximal cliques are the
        right vertices and the left vertices are retained from the given
        graph. Right and left vertices are connected if the bottom vertex
        belongs to the clique represented by a top vertex.

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
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
        from bipartite_graph import BipartiteGraph
        import networkx
        return BipartiteGraph(networkx.make_clique_bipartite(self.networkx_graph(copy=False), **kwds))

    def independent_set(self):
        """
        Returns a maximal independent set, which is a set of vertices which
        induces an empty subgraph. Uses Cliquer [NisOst2003]_.

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
            to convert a digraph to an undirected graph.

        EXAMPLES::

            sage: C=graphs.PetersenGraph()
            sage: C.independent_set()
            [0, 3, 6, 7]
        """
        from sage.graphs.cliquer import max_clique
        return max_clique(self.complement())

    def cliques_vertex_clique_number(self, algorithm="cliquer", vertices=None,
                                     cliques=None):
        """
        Returns a dictionary of sizes of the largest maximal cliques containing
        each vertex, keyed by vertex. (Returns a single value if only one
        input vertex).

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
            to convert a digraph to an undirected graph.

        INPUT:

         - ``algorithm`` - either ``cliquer`` or ``networkx``

           - ``cliquer`` - This wraps the C program Cliquer [NisOst2003]_.

           - ``networkx`` - This function is based on NetworkX's implementation
                of the Bron and Kerbosch Algorithm [BroKer1973]_.

        -  ``vertices`` - the vertices to inspect (default is entire graph).
           Ignored unless ``algorithm=='networkx'``.

        -  ``cliques`` - list of cliques (if already computed).  Ignored unless
           ``algorithm=='networkx'``.

        EXAMPLES::

            sage: C = Graph('DJ{')
            sage: C.cliques_vertex_clique_number()
            {0: 2, 1: 4, 2: 4, 3: 4, 4: 4}
            sage: E = C.cliques_maximal()
            sage: E
            [[4, 0], [4, 1, 2, 3]]
            sage: C.cliques_vertex_clique_number(cliques=E,algorithm="networkx")
            {0: 2, 1: 4, 2: 4, 3: 4, 4: 4}
            sage: F = graphs.Grid2dGraph(2,3)
            sage: X = F.cliques_vertex_clique_number(algorithm="networkx")
            sage: for v in sorted(X.iterkeys()):
            ...    print v, X[v]
            (0, 0) 2
            (0, 1) 2
            (0, 2) 2
            (1, 0) 2
            (1, 1) 2
            (1, 2) 2
            sage: F.cliques_vertex_clique_number(vertices=[(0, 1), (1, 2)])
            {(0, 1): 2, (1, 2): 2}
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_vertex_clique_number()
            {0: 3, 1: 3, 2: 3, 3: 3}

        """

        if algorithm=="cliquer":
            from sage.graphs.cliquer import clique_number
            if vertices==None:
                vertices=self
            value={}
            for v in vertices:
                value[v] = 1+clique_number(self.subgraph(self.neighbors(v)))
                self.subgraph(self.neighbors(v)).plot()
            return value
        elif algorithm=="networkx":
            import networkx
            return networkx.node_clique_number(self.networkx_graph(copy=False),vertices, cliques)
        else:
            raise NotImplementedError("Only 'networkx' and 'cliquer' are supported.")

    def cliques_containing_vertex(self, vertices=None, cliques=None):
        """
        Returns the cliques containing each vertex, represented as a dictionary
        of lists of lists, keyed by vertex. (Returns a single list if only one
        input vertex).

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
            to convert a digraph to an undirected graph.

        INPUT:

        -  ``vertices`` - the vertices to inspect (default is
           entire graph)

        -  ``cliques`` - list of cliques (if already
           computed)

        EXAMPLES::

            sage: C = Graph('DJ{')
            sage: C.cliques_containing_vertex()
            {0: [[4, 0]], 1: [[4, 1, 2, 3]], 2: [[4, 1, 2, 3]], 3: [[4, 1, 2, 3]], 4: [[4, 0], [4, 1, 2, 3]]}
            sage: E = C.cliques_maximal()
            sage: E
            [[4, 0], [4, 1, 2, 3]]
            sage: C.cliques_containing_vertex(cliques=E)
            {0: [[4, 0]], 1: [[4, 1, 2, 3]], 2: [[4, 1, 2, 3]], 3: [[4, 1, 2, 3]], 4: [[4, 0], [4, 1, 2, 3]]}
            sage: F = graphs.Grid2dGraph(2,3)
            sage: X = F.cliques_containing_vertex()
            sage: for v in sorted(X.iterkeys()):
            ...    print v, X[v]
            (0, 0) [[(0, 1), (0, 0)], [(1, 0), (0, 0)]]
            (0, 1) [[(0, 1), (0, 0)], [(0, 1), (0, 2)], [(0, 1), (1, 1)]]
            (0, 2) [[(0, 1), (0, 2)], [(1, 2), (0, 2)]]
            (1, 0) [[(1, 0), (0, 0)], [(1, 0), (1, 1)]]
            (1, 1) [[(0, 1), (1, 1)], [(1, 2), (1, 1)], [(1, 0), (1, 1)]]
            (1, 2) [[(1, 2), (0, 2)], [(1, 2), (1, 1)]]
            sage: F.cliques_containing_vertex(vertices=[(0, 1), (1, 2)])
            {(0, 1): [[(0, 1), (0, 0)], [(0, 1), (0, 2)], [(0, 1), (1, 1)]], (1, 2): [[(1, 2), (0, 2)], [(1, 2), (1, 1)]]}
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_containing_vertex()
            {0: [[0, 1, 2], [0, 1, 3]], 1: [[0, 1, 2], [0, 1, 3]], 2: [[0, 1, 2]], 3: [[0, 1, 3]]}

        """
        import networkx
        return networkx.cliques_containing_node(self.networkx_graph(copy=False),vertices, cliques)

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
            [(3, 4, {}), (2, 4, {}), (1, 4, {}), (0, 4, {})]
            sage: g.min_spanning_tree(algorithm='Prim edge', starting_vertex=2, weight_function=weight)
            [(2, 4, {}), (3, 4, {}), (1, 4, {}), (0, 4, {})]
            sage: g.min_spanning_tree(algorithm='Prim fringe', starting_vertex=2, weight_function=weight)
            [(2, 4), (4, 3), (4, 1), (4, 0)]
        """
        if self.is_connected()==False:
            return False

        if algorithm=='Kruskal':
            # Kruskal's algorithm
            edges=[]
            sorted_edges_iterator=iter(sorted(self.edges(), key=weight_function))
            union_find = dict()
            while len(edges) < self.order()-1:
                # get next edge
                e=sorted_edges_iterator.next()
                components=[]
                for start_v in e[0:2]:
                    v=start_v
                    children=[]

                    # Find the component a vertex lives in.
                    while union_find.has_key(v):
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

            fringe_list=dict([u,(weight_function((v,u)),v)] for u in self[v])

            cmp_fun = lambda x: fringe_list[x]

            for i in xrange(self.order()-1):
                # Find the smallest-weight fringe vertex
                u=min(fringe_list,key=cmp_fun)
                edges.append((fringe_list[u][1],u))
                tree.add(u)
                fringe_list.pop(u)

                # Update fringe list
                for neighbor in [v for v in self[u] if v not in tree]:
                    w=weight_function((u,neighbor))
                    if neighbor not in fringe_list or fringe_list[neighbor][0]>w:
                        fringe_list[neighbor]=(w,u)
            return edges

        elif algorithm=='Prim edge':
            if starting_vertex is None:
                v = self.vertex_iterator().next()
            else:
                v = starting_vertex
            sorted_edges=sorted(self.edges(), key=weight_function)
            tree=set([v])
            edges=[]

            for _ in xrange(self.order()-1):
                # Find a minimum-weight edge connecting a vertex in
                # the tree to something outside the tree.  Remove the
                # edges between tree vertices for efficiency.
                i=0
                while True:
                    e = sorted_edges[i]
                    v0,v1=e[0],e[1]
                    if v0 in tree:
                        del sorted_edges[i]
                        if v1 in tree: continue
                        edges.append(e)
                        tree.add(v1)
                        break
                    elif v1 in tree:
                        del sorted_edges[i]
                        edges.append(e)
                        tree.add(v0)
                        break
                    else:
                        i += 1
            return edges

        elif algorithm=="networkx":
            import networkx
            G = networkx.Graph([(u,v,dict(weight=weight_function((u,v)))) for u,v,l in self.edge_iterator()])
            return list(networkx.mst(G))
        else:
            raise NotImplementedError, "Minimum Spanning Tree algorithm '%s' is not implemented."%algorithm

    def modular_decomposition(self):
        r"""
        Returns the modular decomposition corresponding
        to the current graph.

        Crash course on modular decomposition:

        A module `M` of a graph `G` is a proper subset of its vertices
        such that for all `u \in V(G)-M, v,w\in M` the relation `u
        \sim v \Leftrightarrow u \sim w` holds, where `\sim` denotes
        the adjacency relation in `G`. Equivalently, `M \subset V(G)`
        is a module if all its vertices have the same adjacency
        relations with each vertex outside of the module (vertex by
        vertex).

        Hence, for a set like a module, it is very easy to encode the
        information of the adjacencies between the vertices inside and
        outside the module -- we can actually add a new vertex `v_M`
        to our graph representing our module `M`, and let `v_M` be
        adjacent to `u\in V(G)-M` if and only if some `v\in M` (and
        hence all the vertices contained in the module) is adjacent to
        `u`. We can now independently (and recursively) study the
        structure of our module `M` and the new graph `G-M+\{v_M\}`,
        without any loss of information.

        Here are two very simple modules :

            * A connected component `C` (or the union of some --but
              not all-- of them) of a disconnected graph `G`, for
              instance, is a module, as no vertex of `C` has a
              neighbor outside of it.

            * An anticomponent `C` (or the union of some --but not
              all-- of them) of an non-anticonnected graph `G`, for
              the same reason (it is just the complement of the
              previous graph !).

        These modules being of special interest, the disjoint union of
        graphs is called a Parallel composition, and the complement of
        a disjoint union is called a Parallel composition. A graph
        whose only modules are singletons is called Prime.

        For more information on modular decomposition, in particular
        for an explanation of the terms "Parallel," "Prime" and
        "Serie," see the `Wikipedia article on modular decomposition
        <http://en.wikipedia.org/wiki/Modular_decomposition>`_.

        You may also be interested in the survey from Michel Habib and
        Christophe Paul entitled "A survey on Algorithmic aspects of
        modular decomposition" [HabPau10]_.

        OUTPUT:

        A pair of two values (recursively encoding the decomposition) :

            * The type of the current module :

                * ``"Parallel"``
                * ``"Prime"``
                * ``"Serie"``

            * The list of submodules (as list of pairs ``(type, list)``,
              recursively...) or the vertex's name if the module is a
              singleton.

        EXAMPLES:

        The Bull Graph is prime::

            sage: graphs.BullGraph().modular_decomposition()
            ('Prime', [3, 4, 0, 1, 2])

        The Petersen Graph too::

            sage: graphs.PetersenGraph().modular_decomposition()
            ('Prime', [2, 6, 3, 9, 7, 8, 0, 1, 5, 4])

        This a clique on 5 vertices with 2 pendant edges, though, has a more
        interesting decomposition ::

            sage: g = graphs.CompleteGraph(5)
            sage: g.add_edge(0,5)
            sage: g.add_edge(0,6)
            sage: g.modular_decomposition()
            ('Serie', [0, ('Parallel', [5, ('Serie', [1, 4, 3, 2]), 6])])

        ALGORITHM:

        This function uses a C implementation of a 2-step algorithm
        implemented by Fabien de Montgolfier [FMDec]_ :

            * Computation of a factorizing permutation [HabibViennot1999]_.

            * Computation of the tree itself [CapHabMont02]_.

        .. SEEALSO::

        - :meth:`is_prime` -- Tests whether a graph is prime.

        REFERENCE:

        .. [FMDec] Fabien de Montgolfier
          http://www.liafa.jussieu.fr/~fm/algos/index.html

        .. [HabibViennot1999] Michel Habib, Christiphe Paul, Laurent Viennot
          Partition refinement techniques: An interesting algorithmic tool kit
          International Journal of Foundations of Computer Science
          vol. 10 n2 pp.147--170, 1999

        .. [CapHabMont02] C. Capelle, M. Habib et F. de Montgolfier
          Graph decomposition and Factorising Permutations
          Discrete Mathematics and Theoretical Computer Sciences, vol 5 no. 1 , 2002.

        .. [HabPau10] Michel Habib and Christophe Paul
          A survey of the algorithmic aspects of modular decomposition
          Computer Science Review
          vol 4, number 1, pages 41--59, 2010
          http://www.lirmm.fr/~paul/md-survey.pdf
        """

        from sage.graphs.modular_decomposition.modular_decomposition import modular_decomposition

        D = modular_decomposition(self)

        id_label = dict(enumerate(self.vertices()))

        relabel = lambda x : (x[0], map(relabel,x[1])) if isinstance(x,tuple) else id_label[x]

        return relabel(D)

    def is_prime(self):
        r"""
        Tests whether the current graph is prime. A graph is prime if
        all its modules are trivial (i.e. empty, all of the graph or
        singletons)-- see `self.modular_decomposition?`.

        EXAMPLE:

        The Petersen Graph and the Bull Graph are both prime ::

            sage: graphs.PetersenGraph().is_prime()
            True
            sage: graphs.BullGraph().is_prime()
            True

        Though quite obviously, the disjoint union of them is not::

            sage: (graphs.PetersenGraph() + graphs.BullGraph()).is_prime()
            False
        """

        D = self.modular_decomposition()

        return D[0] == "Prime" and len(D[1]) == self.order()

    def _gomory_hu_tree(self, vertices=None, method="FF"):
        r"""
        Returns a Gomory-Hu tree associated to self.

        This function is the private counterpart of ``gomory_hu_tree()``,
        with the difference that it has an optional argument
        needed for recursive computations, which the user is not
        interested in defining himself.

        See the documentation of ``gomory_hu_tree()`` for more information.

        INPUT:

        - ``vertices`` - a set of "real" vertices, as opposed to the
          fakes one introduced during the computations. This variable is
          useful for the algorithm and for recursion purposes.

        - ``method`` -- There are currently two different
          implementations of this method :

              * If ``method = "FF"`` (default), a Python
                implementation of the Ford-Fulkerson algorithm is
                used.

              * If ``method = "LP"``, the flow problem is solved using
                Linear Programming.

        EXAMPLE:

        This function is actually tested in ``gomory_hu_tree()``, this
        example is only present to have a doctest coverage of 100%.

            sage: g = graphs.PetersenGraph()
            sage: t = g._gomory_hu_tree()
        """
        from sage.sets.set import Set

        # The default capacity of an arc is 1
        from sage.rings.real_mpfr import RR
        capacity = lambda label: label if label in RR else 1

        # Keeping the graph's embedding
        pos = False

        # Small case, not really a problem ;-)
        if self.order() == 1:
            return self.copy()

        # This is a sign that this is the first call
        # to this recursive function
        if vertices is None:
            # Now is the time to care about positions
            pos = self.get_pos()

            # if the graph is not connected, returns the union
            # of the Gomory-Hu tree of each component
            if not self.is_connected():
                g = Graph()
                for cc in self.connected_components_subgraphs():
                    g = g.union(cc._gomory_hu_tree(method=method))
                g.set_pos(self.get_pos())
                return g
            # All the vertices is this graph are the "real ones"
            vertices = Set(self.vertices())

        # There may be many vertices, though only one which is "real"
        if len(vertices) == 1:
            g = Graph()
            g.add_vertex(vertices[0])
            return g

        # Take any two vertices
        u,v = vertices[0:2]

        # Recovers the following values
        # flow is the connectivity between u and v
        # edges of a min cut
        # sets1, sets2 are the two sides of the edge cut
        flow,edges,[set1,set2] = self.edge_cut(u, v, use_edge_labels=True, vertices=True, method=method)

        # One graph for each part of the previous one
        g1,g2 = self.subgraph(set1), self.subgraph(set2)

        # Adding the fake vertex to each part
        g1_v = Set(set2)
        g2_v = Set(set1)
        g1.add_vertex(g1_v)
        g1.add_vertex(g2_v)

        # Each part of the graph had many edges going to the other part
        # Now that we have a new fake vertex in each part
        # we just say that the edges which were in the cut and going
        # to the other side are now going to this fake vertex

        # We must preserve the labels. They sum.

        for e in edges:
            x,y = e[0],e[1]
            # Assumes x is in g1
            if x in g2:
                x,y = y,x
            # If the edge x-g1_v exists, adds to its label the capacity of arc xy
            if g1.has_edge(x, g1_v):
                g1.set_edge_label(x, g1_v, g1.edge_label(x, g1_v) + capacity(self.edge_label(x, y)))
            else:
                # Otherwise, creates it with the good label
                g1.add_edge(x, g1_v, capacity(self.edge_label(x, y)))
            # Same thing for g2
            if g2.has_edge(y, g2_v):
                g2.set_edge_label(y, g2_v, g2.edge_label(y, g2_v) + capacity(self.edge_label(x, y)))
            else:
                g2.add_edge(y, g2_v, capacity(self.edge_label(x, y)))

        # Recursion for the two new graphs... The new "real" vertices are the intersection with
        # with the previous set of "real" vertices
        g1_tree = g1._gomory_hu_tree(vertices=(vertices & Set(g1.vertices())), method=method)
        g2_tree = g2._gomory_hu_tree(vertices=(vertices & Set(g2.vertices())), method=method)

        # Union of the two partial trees ( it is disjoint, but
        # disjoint_union does not preserve the name of the vertices )
        g = g1_tree.union(g2_tree)

        # An edge to connect them, with the appropriate label
        g.add_edge(g1_tree.vertex_iterator().next(), g2_tree.vertex_iterator().next(), flow)

        if pos:
            g.set_pos(pos)

        return g

    def gomory_hu_tree(self, method="FF"):
        r"""
        Returns a Gomory-Hu tree of self.

        Given a tree `T` with labeled edges representing capacities, it is very
        easy to determine the maximal flow between any pair of vertices :
        it is the minimal label on the edges of the unique path between them.

        Given a graph `G`, a Gomory-Hu tree `T` of `G` is a tree
        with the same set of vertices, and such that the maximal flow
        between any two vertices is the same in `G` as in `T`. See the
        `Wikipedia article on Gomory-Hu tree <http://en.wikipedia.org/wiki/Gomory%E2%80%93Hu_tree>`_.
        Note that, in general, a graph admits more than one Gomory-Hu tree.

        INPUT:

        - ``method`` -- There are currently two different
          implementations of this method :

              * If ``method = "FF"`` (default), a Python
                implementation of the Ford-Fulkerson algorithm is
                used.

              * If ``method = "LP"``, the flow problems are solved
                using Linear Programming.

        OUTPUT:

        graph with labeled edges

        EXAMPLE:

        Taking the Petersen graph::

            sage: g = graphs.PetersenGraph()
            sage: t = g.gomory_hu_tree()

        Obviously, this graph is a tree::

            sage: t.is_tree()
            True

        Note that if the original graph is not connected, then the
        Gomory-Hu tree is in fact a forest::

            sage: (2*g).gomory_hu_tree().is_forest()
            True
            sage: (2*g).gomory_hu_tree().is_connected()
            False

        On the other hand, such a tree has lost nothing of the initial
        graph connectedness::

            sage: all([ t.flow(u,v) == g.flow(u,v) for u,v in Subsets( g.vertices(), 2 ) ])
            True

        Just to make sure, we can check that the same is true for two vertices
        in a random graph::

            sage: g = graphs.RandomGNP(20,.3)
            sage: t = g.gomory_hu_tree()
            sage: g.flow(0,1) == t.flow(0,1)
            True

        And also the min cut::

            sage: g.edge_connectivity() == min(t.edge_labels())
            True
        """
        return self._gomory_hu_tree(method=method)


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
            sage: classes = g.two_factor_petersen()
            sage: for c in classes:
            ...     gg = Graph()
            ...     gg.add_edges(c)
            ...     print max(gg.degree())<=2
            True
            True
            True
            sage: Set(set(classes[0]) | set(classes[1]) | set(classes[2])).cardinality() == g.size()
            True

        ::

            sage: g = graphs.CirculantGraph(24, [7, 11])
            sage: cl = g.two_factor_petersen()
            sage: g.plot(edge_colors={'black':cl[0], 'red':cl[1]})

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
        classes = edge_coloring(g)

        # The edges in the classes are of the form ((-1,u),(1,v))
        # and have to be translated back to (u,v)
        classes_b = []
        for c in classes:
            classes_b.append([(u,v) for ((uu,u),(vv,v)) in c])

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

