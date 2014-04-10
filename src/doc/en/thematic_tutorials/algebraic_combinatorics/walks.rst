===============
Walks in graphs
===============

This section provides some examples on Chapter 1 of Stanley's book.

We begin by creating a graph with 4 vertices::

    sage: G = Graph(4)
    sage: G
    Graph on 4 vertices

This graph has no edges yet!
::

    sage: G.vertices()
    [0, 1, 2, 3]
    sage: G.edges()
    []

Before we can add edges, we need to tell Sage that our graph can
have loops and multiple edges.::

    sage: G.allow_loops(True)
    sage: G.allow_multiple_edges(True)

Now we are ready to add our edges by specifying a tuple of vertices that
are connected by an edge. If there are multiple edges, we need to add
the tuple with multiplicity.::

    sage: G.add_edges([(0,0),(0,0),(0,1),(0,3),(1,3),(1,3)])

Now let us look at the graph!

.. image:: ../media/graph.png
   :scale: 75
   :align: center
