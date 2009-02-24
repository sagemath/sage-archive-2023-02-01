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
``C.show(vertex_labels=False, node_size=60, graph_border=True, figsize=[9,8])``
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

As another option to do graph theory in Sage, you may also load
Leonard Soicher's GAP package {GRAPE}
(http://www.gap-system.org/Packages/grape.html), which in turn
calls the C programs in Brendan McKay's ``nauty``
(http://cs.anu.edu.au/people/bdm/nauty/). These packages require a
UNIX environment to be installed (such as Linux or Cygwin - see the
GRAPE readme file http://www.maths.qmul.ac.uk/~leonard/grape/README
file for details). To install GRAPE in Sage, see :ref:`section-installAll`.

.. index::
   pair: graph; Cayley
   single: GRAPE

::

    sage: print gap.eval('LoadPackage("grape")')   # need optional gap packages
    true
    sage: print gap.eval("C := CayleyGraph(SymmetricGroup(4),[(1,2),(2,3),(3,4)])")               # optional gap package
    rec( isGraph := true, order := 24,
      group := Group([ (1,10,17,19)(2,9,18,20)(3,12,14,21)(4,11,13,22)(5,7,16,
            23)(6,8,15,24), (1,7)(2,8)(3,9)(4,10)(5,11)(6,12)(13,15)(14,16)(17,
            18)(19,21)(20,22)(23,24) ]),
      schreierVector := [ -1, 1, 1, 2, 2, 1, 2, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1, 2,
          1, 1, 2, 2, 1, 2 ], adjacencies := [ [ 2, 3, 7 ] ],
      representatives := [ 1 ],
      names := [ (), (3,4), (2,3), (2,3,4), (2,4,3), (2,4), (1,2), (1,2)(3,4),
          (1,2,3), (1,2,3,4), (1,2,4,3), (1,2,4), (1,3,2), (1,3,4,2), (1,3),
          (1,3,4), (1,3)(2,4), (1,3,2,4), (1,4,3,2), (1,4,2), (1,4,3), (1,4),
          (1,4,2,3), (1,4)(2,3) ], isSimple := true )
    sage: print gap.eval("Girth(C)")                   # optional gap package
    4
    sage: print gap.eval("Diameter(C)")             # optional gap package
    6
    sage: print gap.eval("AutGroupGraph(C)")  # optional gap package (uses nauty)
    Group([ (2,7)(4,13)(5,9)(6,15)(10,19)(12,21)(16,20)(18,23),
      (1,2)(3,5)(4,6)(7,8)(9,11)(10,12)(13,19)(14,20)(15,21)(16,22)(17,23)(18,24),
      (1,3)(2,4)(5,6)(7,13)(8,14)(9,15)(10,16)(11,17)(12,18)(19,20)(21,23)(22,24)
     ])

The command ``AutGroupGraph`` uses the 2.2 version of ``nauty``.

Here is a more Pythonic version (thanks to Jack Schmidt for helping
with this):

::

    C3 = CyclicPermutationGroup(3)._gap_()
    C = C3.CayleyGraph(C3.GeneratorsOfGroup())
    V = C.Vertices().Elements()
    E = C.UndirectedEdges()
    L = [[y[1] for y in [x for x in E if v in x] if y[1]!=v]+[y[2]
    for y in [x for x in E if v in x] if y[2]!=v] for v in V]
    d = dict(zip(V,L))
    G = Graph(d)
    show(G.plot())

.. index::
   pair: graph; adjacency matrix

.. section_adjacency:

Graphs from adjacency matrices
==============================

To construct the graph Gamma with :math:`n \times n` adjacency
matrix :math:`A`, you want a graph :math:`X` so that the
vertex-set of Gamma is :math:`\{1,..., n\}`, and :math:`[i,j]`
is an edge of Gamma if and only if :math:`A[i][j] = 1`.

Here is an example of the syntax in (copied from Robert Miller's
SageDays 3 talk):

::

    sage: M = Matrix ([ [0, 1, 1],
    ...   [1, 0, 1],
    ...   [1, 1, 0] ])
    ...   # (the order is the number of edges)
    sage: G = Graph(M); G.order()
    3



Here is an example of the syntax in GRAPE:

::

    sage: print gap.eval("A := [[0,1,0],[1,0,0],[0,0,1]]")
    [ [ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ]
    sage: print gap.eval("G := Group( (1,2) )")
    Group([ (1,2) ])
    sage: print gap.eval("Gamma := Graph( G, [1..3], OnPoints, function(x,y) return A[x][y] = 1; end,true )")  # optional gap package
    rec( isGraph := true, order := 3, group := Group([ (1,2) ]),
      schreierVector := [ -1, 1, -2 ], adjacencies := [ [ 2 ], [ 3 ] ],
      representatives := [ 1, 3 ], names := [ 1, 2, 3 ] )
    sage: print gap.eval("Adjacency(Gamma,1)")           # optional gap package
    [ 2 ]
    sage: print gap.eval("Adjacency(Gamma,2)")           # optional gap package
    [ 1 ]
    sage: print gap.eval("Adjacency(Gamma,3)")           # optional gap package
    [ 3 ]
    sage: print gap.eval("IsEdge( Gamma, [ 1, 2 ] )")    # optional gap package
    true
    sage: print gap.eval("Vertices( Gamma )")            # optional gap package
    [ 1 .. 3 ]

Define the distance :math:`d(x,y)` from :math:`x` to
:math:`y` to be the minimum length of a (directed) path in Gamma
joining a vertex :math:`x` to a vertex :math:`y` if such a path
exists, and :math:`-1` otherwise.

::

    sage: print gap.eval("Distance( Gamma, 1, 3 )")      # optional gap package
    -1

A diameter of :math:`-1` is returned if Gamma is not (strongly)
connected. Otherwise, the diameter of Gamma is equal to the maximum
(directed) distance :math:`d(x,y)` in gamma (as :math:`x` and
:math:`y` range over all the vertices of Gamma).

::

    sage: print gap.eval("Distance( Gamma, 1, 2 )")      # optional gap package
    1
