Graph Theory in Sage
====================

First steps
-----------

A Graph is a set of vertices connected by edges
(see `Wikipedia <http://en.wikipedia.org/wiki/Graph_(mathematics)>`_).

One can very easily create a graph in sage by typing::

    sage: g=Graph()

By typing the name of the Graph, one can get some basic information
about it::

    sage: g
    Graph on 0 vertices

Sage contains a large collection of predefined graph classes
that can be listed this way:

#. type in Sage: ``graphs.`` (do not press "Enter", and do not
   forget the final "."), then

#. hit the tab key two times in a row

You will see the list of methods defined in the class ``graphs``,
all of which generate graphs you can play with!

If you want to see what they look like, you can plot the graphs::

    sage: g=graphs.PetersenGraph()
    sage: g.plot()
    sage: g=graphs.ChvatalGraph()
    sage: g.plot()

If you are curious about what these graphs are, for example, if you
wonder what ``RandomGNP`` actually is, use the standard Sage help::

    sage: graphs.RandomGNP?

Once you have defined the graph you want, you can begin
to work on it by using the hundreds of functions for graphs
in the Sage library.  If your graph is named ``g``, you can
list these functions using tab completion:

#. type in Sage: ``g.`` (do not press "Enter", and do not forget
   the final "."), then

#. hit the tab key two times in a row

As usual, you can get some information about what these
functions do by typing (if you want to know about the :meth:`diameter()`
method)::

    sage: g.diameter?

To find out if a graph is connected, use :meth:`is_connected`::

   sage: g=graphs.RandomGNP(100, 0.01)
   sage: g.is_connected()

You can plot each connected component separately::

    sage: for component in g.connected_components():
    ...      g.subgraph(component).plot()


Basic notions
-------------

The first elements you want to find in a graph are its vertices
and its edges. The vertices of a graph ``g`` are returned
by :meth:`g.vertices` and its edges are returned by :meth:`g.edges`.
In Sage, the edges of a graph are represented as tuples ``(u,v,l)``
where ``u`` and ``v`` are vertices, and ``l`` a label attached
to the edge (this label can be of any type, even though
several functions expect it to be a real value).

The methods :meth:`g.order` and :meth:`g.size` respectively return the number
of vertices and the number of edges.

At any moment, you can display the adjacency matrix of your graph
by using the method :meth:`g.adjacency_matrix` and plot the graph with
:meth:`g.plot`.

What interesting things can you do with Graphs in Sage?
-------------------------------------------------------

Compute maximum matchings
^^^^^^^^^^^^^^^^^^^^^^^^^

Maximum Matchings, a polynomial problem in Graph Theory, has many
different applications.  In the subsections that follow, we will look
at a few of these applications.

For more information on matchings, see `Matching
<http://en.wikipedia.org/wiki/Matching>`_.

Small company
"""""""""""""

Let us say that you are in charge of a small company with 4 employees
`\{e_1,e_2,e_3,e_4\}` and have several tasks `\{t_1,t_2,t_3,t_4\}`
to give them. Unfortunately, no worker is skilled enough to do all of them:

* `e_1` can do `t_1, t_3, t_4`
* `e_2` can do `t_1, t_3, t_5`
* `e_3` can do `t_1, t_2, t_3, t_4, t_5`
* `e_4` can do `t_4, t_5`
* `e_5` can do `t_2, t_4`

You are lucky if you do not know how to solve this problem manually,
because this is typically an application of matching in graphs
(and if you have found the solution, I assure you it gets harder
when you have more employees).

To solve this problem, create the graph corresponding to the
information above, and solve the matching problem::

    sage: g=Graph({"e0":['t1', 't3', 't4'],"e1":['t1', 't3',
            't5'],"e2":['t1', 't2', 't3', 't4', 't5'],
            "e3":['t4', 't5'],"e4":['t2', 't4']})
    sage: print g.max_matching()
    [('e2', 't4', None), ('e3', 't5', None), ('e0', 't3', None),
     ('e1', 't1', None), ('e4', 't2', None)]

If you prefer to "see" the result, you can also type::

    sage: g.plot(edge_colors={"red":g.max_matching()})

Wasn't that simple?

Summer camp
"""""""""""

You now have under your responsibility five rooms and ten children
`\{c_0, \dots ,c_9\}`. You need to decide which of them will
sleep in the same rooms, but you do not want two of them to be
together if they do not like each other or if you expect trouble
from the pair. Here are the constraints:

* `c_0` can sleep with `c_5`
* `c_1` can sleep with `c_5, c_8`
* `c_2` can sleep with `c_3, c_8, c_9`
* `c_3` can sleep with `c_9`
* `c_4` can sleep with `c_9`
* `c_5` can sleep with `c_9`
* `c_6` can sleep with `c_7, c_9`
* `c_7` can sleep with `c_9`

As done previously, this defines a graph. Now create it in Sage, and
ask for a maximum matching::

    sage: g=Graph({'c0':['c5'],'c1':['c5', 'c8'],'c2':['c3',
            'c8', 'c9'],'c3':['c9'],'c4':['c9'],'c5':['c9'],
            'c6':['c7', 'c9'],'c7':['c9']})
    sage: print g.max_matching()
    [('c0', 'c5', None), ('c6', 'c7', None), ('c2', 'c3', None),
     ('c4', 'c9', None), ('c1', 'c8', None)]

If you prefer to "see" the result, you can also type::

    sage: g.plot(edge_colors={"red":g.max_matching()})

And this is another problem Sage solved for you!

Vertex coloring
^^^^^^^^^^^^^^^

You are in front of a map of Western Europe that you would like
to color. Obviously, you can not color both France and Italy
with the same color, as they have a common boundary, and you would
not like to mix the two. Actually, you want to color:

* Austria
* Belgium
* France
* Germany
* Ireland
* Italy
* Luxembourg
* Netherlands
* Portugal
* Spain
* Swiss
* United Kingdom

And would like to know how many colors you need, and how to color
them. Well, as Sage was especially built to help you solve this
kind of tremendously exciting question, here is the way to solve them:

#. Create the graph of Western Europe in Sage.
#. Use the :meth:`vertex_coloring` method.

In Sage::

    sage: g=Graph({"France":["Italy","Spain","Swiss","Luxembourg","Belgium",
                             "Germany","Austria"],
                   "Spain":["Portugal"],
                   "Italy":["Swiss","Austria"],
                   "Swiss":["Germany"],
                   "Germany":["Luxembourg","Belgium","Netherlands"],
                   "Belgium":["Luxembourg","Netherlands"],
                   "United Kingdom":["Ireland"]})
    sage: g.vertex_coloring()
    [['France', 'Portugal', 'Netherlands', 'Ireland'],
     ['Germany', 'Spain', 'Austria', 'United Kingdom'],
     ['Belgium', 'Swiss'],
     ['Luxembourg', 'Italy']]

You can now look for your pens---four of them.

For more information on graph
coloring, see `Graph coloring <http://en.wikipedia.org/wiki/Graph_coloring>`_.

For more information on why it could not have required more pens, see the
`Four color theorem
<http://en.wikipedia.org/wiki/Four_color_theorem>`_.

Edge coloring
^^^^^^^^^^^^^

You are organizing a soccer tournament, with ten different teams that
are to play against each other. The teams will play every Wednesday
and will not be able to play two times on the same day. How can you
schedule them in such a way that the tournament will not last for too
long?

This is an easy application of the Edge Coloring problem on a
complete graph. If you number your teams as `1, \dots,10`, here
is how you can obtain your scheduling::

    sage: g=graphs.CompleteGraph(10)
    sage: g.edge_coloring()
    [[(2, 9, None), (3, 7, None), (5, 6, None), (0, 8, None), (1, 4, None)],
     [(3, 5, None), (1, 2, None), (7, 9, None), (0, 6, None), (4, 8, None)],
     [(5, 7, None), (0, 3, None), (1, 6, None), (4, 9, None), (2, 8, None)],
     [(1, 7, None), (0, 9, None), (4, 5, None), (2, 3, None), (6, 8, None)],
     [(2, 6, None), (0, 1, None), (4, 7, None), (5, 8, None), (3, 9, None)],
     [(7, 8, None), (3, 4, None), (1, 5, None), (6, 9, None), (0, 2, None)],
     [(5, 9, None), (0, 4, None), (1, 8, None), (3, 6, None), (2, 7, None)],
     [(0, 5, None), (1, 3, None), (6, 7, None), (2, 4, None), (8, 9, None)],
     [(3, 8, None), (4, 6, None), (1, 9, None), (0, 7, None), (2, 5, None)]]

Each line you see is the set of games being played on a particular
day. If you prefer to plot the result, try::

    sage: g.plot(edge_colors=g.edge_coloring(hex_colors=True))

Pretty, isn't it? Each day has its own color.

Two links for more information:

* `About edge coloring <http://en.wikipedia.org/wiki/Edge_coloring>`_
* `About the scheduling of tournaments <http://en.wikipedia.org/wiki/Round-robin>`_
