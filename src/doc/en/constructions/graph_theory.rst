************
Graph theory
************

See the Sage wiki page http://wiki.sagemath.org/graph_survey for an excellent survey
of exisiting graph theory software.

Networkx
========

Networkx (http://networkx.lanl.gov)
"is a Python package for the creation, manipulation, and study of the structure, dynamics, and functions of complex networks".
More details can also be found on
http://wiki.sagemath.org/graph_survey or in Robert Miller's
SageDays 3 talk.

::

    sage: C = graphs.CubeGraph(4)

Now type
``C.show(vertex_labels=False, vertex_size=60, graph_border=True, figsize=[9,8])``
to view this with some of the options.

The digraph below is a :math:`3`-cycle with vertices
:math:`\{0,1,2\}` and edges :math:`0\rightarrow 1`,
:math:`1\rightarrow 2`, :math:`2\rightarrow 0`:

::

    sage: D = DiGraph( { 0: [1], 1: [2], 2: [0]} )

Type ``D.show()`` to view this.

.. _section-cayley:

Cayley graphs
=============

includes wrappers to many NetworkX commands, written mainly by
Emily Kirkman and Robert Miller. The implementation of Cayley
graphs was written by Bobby Moretti and Robert Miller.

::

    sage: G = DihedralGroup(5)
    sage: C = G.cayley_graph(); C
    Digraph on 10 vertices
    sage: C.diameter()
    3
    sage: C.girth()
    2
    sage: C.automorphism_group().order()
    10
    sage: len(C.edges())
    20


.. index::
   pair: graph; adjacency matrix

.. section_adjacency:

Graphs from adjacency matrices
==============================

To construct the graph G with :math:`n \times n` adjacency
matrix :math:`A`, you want a graph :math:`X` so that the
vertex-set of G is :math:`\{1,..., n\}`, and :math:`[i,j]`
is an edge of G if and only if :math:`A[i][j] = 1`.

Here is an example of the syntax in (copied from Robert Miller's
SageDays 3 talk): Define the distance :math:`d(x,y)` from :math:`x` to
:math:`y` to be the minimum length of a (directed) path in Gamma
joining a vertex :math:`x` to a vertex :math:`y` if such a path
exists, and :math:`-1` otherwise.
A diameter of :math:`-1` is returned if G is not (strongly)
connected. Otherwise, the diameter of G is equal to the maximum
(directed) distance :math:`d(x,y)` in G (as :math:`x` and
:math:`y` range over all the vertices of G).

::

    sage: M = Matrix ([ [0, 1, 1], [1, 0, 1], [1, 1, 0] ])
    sage:  # (the order is the number of edges)
    sage: G = Graph(M); G.order()
    3
    sage: G.distance(0,2)
    1
    sage: G.diameter()
    1

