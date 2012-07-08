# -*- coding: utf-8 -*-

r"""
Common graphs

Usage
=====

To see a list of all graph constructors, type "graphs." and then
press the tab key. The documentation for each constructor includes
information about each graph, which provides a useful reference.


Plotting
========

All graphs (i.e., networks) have an associated Sage
graphics object, which you can display::

    sage: G = graphs.WheelGraph(15)
    sage: P = G.plot()
    sage: P.show() # long time

If you create a graph in Sage using the ``Graph``
command, then plot that graph, the positioning of nodes is
determined using the spring-layout algorithm. For the special graph
constructors, which you get using ``graphs.[tab]``, the
positions are preset. For example, consider the Petersen graph with
default node positioning vs. the Petersen graph constructed by this
database::

    sage: petersen_spring = Graph({0:[1,4,5], 1:[0,2,6], 2:[1,3,7], 3:[2,4,8], 4:[0,3,9], 5:[0,7,8], 6:[1,8,9], 7:[2,5,9], 8:[3,5,6], 9:[4,6,7]})
    sage: petersen_spring.show() # long time
    sage: petersen_database = graphs.PetersenGraph()
    sage: petersen_database.show() # long time

For all the constructors in this database (except the octahedral,
dodecahedral, random and empty graphs), the position dictionary is
filled in, instead of using the spring-layout algorithm.

For further visual examples and explanation, see the docstrings
below, particularly for
:meth:`CycleGraph <GraphGenerators.CycleGraph>`,
:meth:`StarGraph <GraphGenerators.StarGraph>`,
:meth:`WheelGraph <GraphGenerators.WheelGraph>`,
:meth:`CompleteGraph <GraphGenerators.CompleteGraph>`, and
:meth:`CompleteBipartiteGraph <GraphGenerators.CompleteBipartiteGraph>`.


.. _organization:

Organization
============

The constructors available in this database are
organized as follows.


Basic structures
----------------

- :meth:`BarbellGraph <GraphGenerators.BarbellGraph>`
- :meth:`BuckyBall <GraphGenerators.BuckyBall>`
- :meth:`BullGraph <GraphGenerators.BullGraph>`
- :meth:`ButterflyGraph <GraphGenerators.ButterflyGraph>`
- :meth:`CircularLadderGraph <GraphGenerators.CircularLadderGraph>`
- :meth:`ClawGraph <GraphGenerators.ClawGraph>`
- :meth:`CycleGraph <GraphGenerators.CycleGraph>`
- :meth:`DiamondGraph <GraphGenerators.DiamondGraph>`
- :meth:`EmptyGraph <GraphGenerators.EmptyGraph>`
- :meth:`Grid2dGraph <GraphGenerators.Grid2dGraph>`
- :meth:`GridGraph <GraphGenerators.GridGraph>`
- :meth:`HouseGraph <GraphGenerators.HouseGraph>`
- :meth:`HouseXGraph <GraphGenerators.HouseXGraph>`
- :meth:`KrackhardtKiteGraph <GraphGenerators.KrackhardtKiteGraph>`
- :meth:`LadderGraph <GraphGenerators.LadderGraph>`
- :meth:`LCFGraph <GraphGenerators.LCFGraph>`
- :meth:`LollipopGraph <GraphGenerators.LollipopGraph>`
- :meth:`PathGraph <GraphGenerators.PathGraph>`
- :meth:`StarGraph <GraphGenerators.StarGraph>`
- :meth:`ToroidalGrid2dGraph <GraphGenerators.ToroidalGrid2dGraph>`
- :meth:`WheelGraph <GraphGenerators.WheelGraph>`


Platonic solids
---------------

- :meth:`DodecahedralGraph <GraphGenerators.DodecahedralGraph>`
- :meth:`HexahedralGraph <GraphGenerators.HexahedralGraph>`
- :meth:`IcosahedralGraph <GraphGenerators.IcosahedralGraph>`
- :meth:`OctahedralGraph <GraphGenerators.OctahedralGraph>`
- :meth:`TetrahedralGraph <GraphGenerators.TetrahedralGraph>`


Named Graphs
------------

- :meth:`Balaban10Cage <GraphGenerators.Balaban10Cage>`
- :meth:`Balaban11Cage <GraphGenerators.Balaban11Cage>`
- :meth:`BidiakisCube <GraphGenerators.BidiakisCube>`
- :meth:`BiggsSmithGraph <GraphGenerators.BiggsSmithGraph>`
- :meth:`BrinkmannGraph <GraphGenerators.BrinkmannGraph>`
- :meth:`ChvatalGraph <GraphGenerators.ChvatalGraph>`
- :meth:`ClebschGraph <GraphGenerators.ClebschGraph>`
- :meth:`CoxeterGraph <GraphGenerators.CoxeterGraph>`
- :meth:`DoubleStarSnark <GraphGenerators.DoubleStarSnark>`
- :meth:`DesarguesGraph <GraphGenerators.DesarguesGraph>`
- :meth:`DurerGraph <GraphGenerators.DurerGraph>`
- :meth:`DyckGraph <GraphGenerators.DyckGraph>`
- :meth:`EllinghamHorton54Graph <GraphGenerators.EllinghamHorton54Graph>`
- :meth:`EllinghamHorton78Graph <GraphGenerators.EllinghamHorton78Graph>`
- :meth:`ErreraGraph <GraphGenerators.ErreraGraph>`
- :meth:`FlowerSnark <GraphGenerators.FlowerSnark>`
- :meth:`FosterGraph <GraphGenerators.FosterGraph>`
- :meth:`FranklinGraph <GraphGenerators.FranklinGraph>`
- :meth:`FruchtGraph <GraphGenerators.FruchtGraph>`
- :meth:`GoldnerHararyGraph <GraphGenerators.GoldnerHararyGraph>`
- :meth:`GrayGraph <GraphGenerators.GrayGraph>`
- :meth:`GrotzschGraph <GraphGenerators.GrotzschGraph>`
- :meth:`HallJankoGraph <GraphGenerators.HallJankoGraph>`
- :meth:`HararyGraph <GraphGenerators.HararyGraph>`
- :meth:`HarriesGraph <GraphGenerators.HarriesGraph>`
- :meth:`HarriesWongGraph <GraphGenerators.HarriesWongGraph>`
- :meth:`HeawoodGraph <GraphGenerators.HeawoodGraph>`
- :meth:`HerschelGraph <GraphGenerators.HerschelGraph>`
- :meth:`HigmanSimsGraph <GraphGenerators.HigmanSimsGraph>`
- :meth:`HoffmanSingletonGraph <GraphGenerators.HoffmanSingletonGraph>`
- :meth:`HoffmanGraph <GraphGenerators.HoffmanGraph>`
- :meth:`LjubljanaGraph <GraphGenerators.LjubljanaGraph>`
- :meth:`McGeeGraph <GraphGenerators.McGeeGraph>`
- :meth:`MoebiusKantorGraph <GraphGenerators.MoebiusKantorGraph>`
- :meth:`MoserSpindle <GraphGenerators.MoserSpindle>`
- :meth:`PappusGraph <GraphGenerators.PappusGraph>`
- :meth:`PetersenGraph <GraphGenerators.PetersenGraph>`
- :meth:`ShrikhandeGraph <GraphGenerators.ShrikhandeGraph>`
- :meth:`ThomsenGraph <GraphGenerators.ThomsenGraph>`
- :meth:`Tutte12Cage <GraphGenerators.Tutte12Cage>`
- :meth:`TutteCoxeterGraph <GraphGenerators.TutteCoxeterGraph>`
- :meth:`WagnerGraph <GraphGenerators.WagnerGraph>`


Families of graphs
------------------

- :meth:`BalancedTree <GraphGenerators.BalancedTree>`
- :meth:`BubbleSortGraph <GraphGenerators.BubbleSortGraph>`
- :meth:`CirculantGraph <GraphGenerators.CirculantGraph>`
- :meth:`CompleteBipartiteGraph <GraphGenerators.CompleteBipartiteGraph>`
- :meth:`CompleteGraph <GraphGenerators.CompleteGraph>`
- :meth:`CubeGraph <GraphGenerators.CubeGraph>`
- :meth:`FibonacciTree <GraphGenerators.FibonacciTree>`
- :meth:`FriendshipGraph <GraphGenerators.FriendshipGraph>`
- :meth:`FuzzyBallGraph <GraphGenerators.FuzzyBallGraph>`
- :meth:`GeneralizedPetersenGraph <GraphGenerators.GeneralizedPetersenGraph>`
- :meth:`HanoiTowerGraph <GraphGenerators.HanoiTowerGraph>`

- :meth:`HyperStarGraph <GraphGenerators.HyperStarGraph>`
- :meth:`KneserGraph <GraphGenerators.KneserGraph>`
- :meth:`LCFGraph <GraphGenerators.LCFGraph>`
- :meth:`MycielskiGraph <GraphGenerators.MycielskiGraph>`
- :meth:`MycielskiStep <GraphGenerators.MycielskiStep>`
- :meth:`NKStarGraph <GraphGenerators.NKStarGraph>`
- :meth:`NStarGraph <GraphGenerators.NStarGraph>`
- :meth:`OddGraph <GraphGenerators.OddGraph>`
- :meth:`PaleyGraph <GraphGenerators.PaleyGraph>`
- :meth:`line_graph_forbidden_subgraphs <GraphGenerators.line_graph_forbidden_subgraphs>`
- :meth:`PermutationGraph <GraphGenerators.PermutationGraph>`
- :meth:`trees <GraphGenerators.trees>`


Pseudofractal graphs
--------------------

- :meth:`DorogovtsevGoltsevMendesGraph <GraphGenerators.DorogovtsevGoltsevMendesGraph>`


Random graphs
-------------

- :meth:`RandomBarabasiAlbert <GraphGenerators.RandomBarabasiAlbert>`
- :meth:`RandomBipartite <GraphGenerators.RandomBipartite>`
- :meth:`RandomGNM <GraphGenerators.RandomGNM>`
- :meth:`RandomGNP <GraphGenerators.RandomGNP>`
- :meth:`RandomHolmeKim <GraphGenerators.RandomHolmeKim>`
- :meth:`RandomInterval <GraphGenerators.RandomInterval>`
- :meth:`RandomLobster <GraphGenerators.RandomLobster>`
- :meth:`RandomNewmanWattsStrogatz <GraphGenerators.RandomNewmanWattsStrogatz>`
- :meth:`RandomRegular <GraphGenerators.RandomRegular>`
- :meth:`RandomShell <GraphGenerators.RandomShell>`
- :meth:`RandomTree <GraphGenerators.RandomTree>`
- :meth:`RandomTreePowerlaw <GraphGenerators.RandomTreePowerlaw>`


Graphs with a given degree sequence
-----------------------------------

- :meth:`DegreeSequence <GraphGenerators.DegreeSequence>`
- :meth:`DegreeSequenceBipartite <GraphGenerators.DegreeSequenceBipartite>`
- :meth:`DegreeSequenceConfigurationModel <GraphGenerators.DegreeSequenceConfigurationModel>`
- :meth:`DegreeSequenceExpected <GraphGenerators.DegreeSequenceExpected>`
- :meth:`DegreeSequenceTree <GraphGenerators.DegreeSequenceTree>`


Miscellaneous
-------------

- :meth:`WorldMap <GraphGenerators.WorldMap>`


AUTHORS:

- Robert Miller (2006-11-05): initial version, empty, random, petersen

- Emily Kirkman (2006-11-12): basic structures, node positioning for
  all constructors

- Emily Kirkman (2006-11-19): docstrings, examples

- William Stein (2006-12-05): Editing.

- Robert Miller (2007-01-16): Cube generation and plotting

- Emily Kirkman (2007-01-16): more basic structures, docstrings

- Emily Kirkman (2007-02-14): added more named graphs

- Robert Miller (2007-06-08-11): Platonic solids, random graphs,
  graphs with a given degree sequence, random directed graphs

- Robert Miller (2007-10-24): Isomorph free exhaustive generation

- Nathann Cohen (2009-08-12): WorldMap

- Michael Yurko (2009-9-01): added hyperstar, (n,k)-star, n-star, and
  bubblesort graphs

- Anders Jonsson (2009-10-15): added generalized Petersen graphs

- Harald Schilly and Yann Laigle-Chapuy (2010-03-24): added Fibonacci Tree

- Jason Grout (2010-06-04): cospectral_graphs

- Edward Scheinerman (2010-08-11): RandomTree

- Ed Scheinerman (2010-08-21): added Grotzsch graph and Mycielski graphs

- Minh Van Nguyen (2010-11-26): added more named graphs

- Keshav Kini (2011-02-16): added Shrikhande and Dyck graphs

- David Coudert (2012-02-10): new RandomGNP generator
"""

###########################################################################

#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#           Copyright (C) 2009 Michael C. Yurko <myurko@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
###########################################################################

# import from Python standard library
from math import sin, cos, pi

# import from Sage library
import graph
from sage.misc.randstate import current_randstate

class GraphGenerators():
    r"""
    A class consisting of constructors for several common graphs, as
    well as orderly generation of isomorphism class representatives. See the
    section :ref:`organization` for a list of supported constructors.

    A list of all graphs and graph structures (other than isomorphism class
    representatives) in this database is available via tab completion. Type
    "graphs." and then hit the tab key to see which graphs are available.

    The docstrings include educational information about each named
    graph with the hopes that this class can be used as a reference.

    For all the constructors in this class (except the octahedral,
    dodecahedral, random and empty graphs), the position dictionary is
    filled to override the spring-layout algorithm.


    ORDERLY GENERATION::

        graphs(vertices, property=lambda x: True, augment='edges', size=None)

    This syntax accesses the generator of isomorphism class
    representatives. Iterates over distinct, exhaustive
    representatives.

    Also: see the use of the optional nauty package for generating graphs
    at the :meth:`nauty_geng` method.

    INPUT:

    - ``vertices`` -- natural number.

    - ``property`` -- (default: ``lambda x: True``) any property to be
      tested on graphs before generation, but note that in general the
      graphs produced are not the same as those produced by using the
      property function to filter a list of graphs produced by using
      the ``lambda x: True`` default. The generation process assumes
      the property has certain characteristics set by the ``augment``
      argument, and only in the case of inherited properties such that
      all subgraphs of the relevant kind (for ``augment='edges'`` or
      ``augment='vertices'``) of a graph with the property also
      possess the property will there be no missing graphs.  (The
      ``property`` argument is ignored if ``degree_sequence`` is
      specified.)

    - ``augment`` -- (default: ``'edges'``) possible values:

      - ``'edges'`` -- augments a fixed number of vertices by
        adding one edge. In this case, all graphs on exactly ``n=vertices`` are
        generated. If for any graph G satisfying the property, every
        subgraph, obtained from G by deleting one edge but not the vertices
        incident to that edge, satisfies the property, then this will
        generate all graphs with that property. If this does not hold, then
        all the graphs generated will satisfy the property, but there will
        be some missing.

      - ``'vertices'`` -- augments by adding a vertex and
        edges incident to that vertex. In this case, all graphs up to
        ``n=vertices`` are generated. If for any graph G satisfying the
        property, every subgraph, obtained from G by deleting one vertex
        and only edges incident to that vertex, satisfies the property,
        then this will generate all graphs with that property. If this does
        not hold, then all the graphs generated will satisfy the property,
        but there will be some missing.

    - ``size`` -- (default: ``None``) the size of the graph to be generated.

    - ``degree_sequence`` -- (default: ``None``) a sequence of non-negative integers,
      or ``None``. If specified, the generated graphs will have these
      integers for degrees. In this case, property and size are both
      ignored.

    - ``loops`` -- (default: ``False``) whether to allow loops in the graph
      or not.

    - ``implementation`` -- (default: ``'c_graph'``) which underlying
      implementation to use (see ``Graph?``).

    - ``sparse`` -- (default: ``True``) ignored if implementation is not
      ``'c_graph'``.

    - ``copy`` (boolean) -- If set to ``True`` (default)
      this method makes copies of the graphs before returning
      them. If set to ``False`` the method returns the graph it
      is working on. The second alternative is faster, but modifying
      any of the graph instances returned by the method may break
      the function's behaviour, as it is using these graphs to
      compute the next ones : only use ``copy_graph = False`` when
      you stick to *reading* the graphs returned.

    EXAMPLES:

    Print graphs on 3 or less vertices::

        sage: for G in graphs(3, augment='vertices'):
        ...    print G
        Graph on 0 vertices
        Graph on 1 vertex
        Graph on 2 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 2 vertices
        Graph on 3 vertices

    Note that we can also get graphs with underlying Cython implementation::

        sage: for G in graphs(3, augment='vertices', implementation='c_graph'):
        ...    print G
        Graph on 0 vertices
        Graph on 1 vertex
        Graph on 2 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 2 vertices
        Graph on 3 vertices

    Print graphs on 3 vertices.

    ::

        sage: for G in graphs(3):
        ...    print G
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices

    Generate all graphs with 5 vertices and 4 edges.

    ::

        sage: L = graphs(5, size=4)
        sage: len(list(L))
        6

    Generate all graphs with 5 vertices and up to 4 edges.

    ::

        sage: L = list(graphs(5, lambda G: G.size() <= 4))
        sage: len(L)
        14
        sage: graphs_list.show_graphs(L) # long time

    Generate all graphs with up to 5 vertices and up to 4 edges.

    ::

        sage: L = list(graphs(5, lambda G: G.size() <= 4, augment='vertices'))
        sage: len(L)
        31
        sage: graphs_list.show_graphs(L)              # long time

    Generate all graphs with degree at most 2, up to 6 vertices.

    ::

        sage: property = lambda G: ( max([G.degree(v) for v in G] + [0]) <= 2 )
        sage: L = list(graphs(6, property, augment='vertices'))
        sage: len(L)
        45

    Generate all bipartite graphs on up to 7 vertices: (see
    http://oeis.org/classic/A033995)

    ::

        sage: L = list( graphs(7, lambda G: G.is_bipartite(), augment='vertices') )
        sage: [len([g for g in L if g.order() == i]) for i in [1..7]]
        [1, 2, 3, 7, 13, 35, 88]

    Generate all bipartite graphs on exactly 7 vertices::

        sage: L = list( graphs(7, lambda G: G.is_bipartite()) )
        sage: len(L)
        88

    Generate all bipartite graphs on exactly 8 vertices::

        sage: L = list( graphs(8, lambda G: G.is_bipartite()) ) # long time
        sage: len(L)                                            # long time
        303

    Remember that the property argument does not behave as a filter,
    except for appropriately inheritable properties::

        sage: property = lambda G: G.is_vertex_transitive()
        sage: len(list(graphs(4, property)))
        1
        sage: len(filter(property, graphs(4)))
        4
        sage: property = lambda G: G.is_bipartite()
        sage: len(list(graphs(4, property)))
        7
        sage: len(filter(property, graphs(4)))
        7

    Generate graphs on the fly: (see http://oeis.org/classic/A000088)

    ::

        sage: for i in range(0, 7):
        ...    print len(list(graphs(i)))
        1
        1
        2
        4
        11
        34
        156

    Generate all simple graphs, allowing loops: (see
    http://oeis.org/classic/A000666)

    ::

        sage: L = list(graphs(5,augment='vertices',loops=True))               # long time
        sage: for i in [0..5]: print i, len([g for g in L if g.order() == i]) # long time
        0 1
        1 2
        2 6
        3 20
        4 90
        5 544

    Generate all graphs with a specified degree sequence (see
    http://oeis.org/classic/A002851)::

        sage: for i in [4,6,8]:  # long time (4s on sage.math, 2012)
        ...       print i, len([g for g in graphs(i, degree_sequence=[3]*i) if g.is_connected()])
        4 1
        6 2
        8 5
        sage: for i in [4,6,8]:  # long time (7s on sage.math, 2012)
        ...       print i, len([g for g in graphs(i, augment='vertices', degree_sequence=[3]*i) if g.is_connected()])
        4 1
        6 2
        8 5

    ::

        sage: print 10, len([g for g in graphs(10,degree_sequence=[3]*10) if g.is_connected()]) # not tested
        10 19

    Make sure that the graphs are really independent and the generator
    survives repeated vertex removal (trac 8458)::

        sage: for G in graphs(3):
        ...       G.delete_vertex(0)
        ...       print(G.order())
        2
        2
        2
        2

    REFERENCE:

    - Brendan D. McKay, Isomorph-Free Exhaustive generation.  *Journal
      of Algorithms*, Volume 26, Issue 2, February 1998, pages 306-324.
    """

    #######################################################################
    #   Basic Structures
    #######################################################################

    def BarbellGraph(self, n1, n2):
        r"""
        Returns a barbell graph with ``2*n1 + n2`` nodes. The argument ``n1``
        must be greater than or equal to 2.

        A barbell graph is a basic structure that consists of a path graph
        of order ``n2`` connecting two complete graphs of order ``n1`` each.

        This constructor depends on `NetworkX <http://networkx.lanl.gov>`_
        numeric labels. In this case, the ``n1``-th node connects to the
        path graph from one complete graph and the ``n1 + n2 + 1``-th node
        connects to the path graph from the other complete graph.

        INPUT:

        - ``n1`` -- integer `\geq 2`. The order of each of the two
          complete graphs.

        - ``n2`` -- nonnegative integer. The order of the path graph
          connecting the two complete graphs.

        OUTPUT:

        A barbell graph of order ``2*n1 + n2``. A ``ValueError`` is
        returned if ``n1 < 2`` or ``n2 < 0``.

        ALGORITHM:

        Uses `NetworkX <http://networkx.lanl.gov>`_.

        PLOTTING:

        Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, each barbell
        graph will be displayed with the two complete graphs in the
        lower-left and upper-right corners, with the path graph connecting
        diagonally between the two. Thus the ``n1``-th node will be drawn at a
        45 degree angle from the horizontal right center of the first
        complete graph, and the ``n1 + n2 + 1``-th node will be drawn 45
        degrees below the left horizontal center of the second complete graph.

        EXAMPLES:

        Construct and show a barbell graph ``Bar = 4``, ``Bells = 9``::

            sage: g = graphs.BarbellGraph(9, 4); g
            Barbell graph: Graph on 22 vertices
            sage: g.show() # long time

        An ``n1 >= 2``, ``n2 >= 0`` barbell graph has order ``2*n1 + n2``. It
        has the complete graph on ``n1`` vertices as a subgraph. It also has
        the path graph on ``n2`` vertices as a subgraph. ::

            sage: n1 = randint(2, 2*10^2)
            sage: n2 = randint(0, 2*10^2)
            sage: g = graphs.BarbellGraph(n1, n2)
            sage: v = 2*n1 + n2
            sage: g.order() == v
            True
            sage: K_n1 = graphs.CompleteGraph(n1)
            sage: P_n2 = graphs.PathGraph(n2)
            sage: s_K = g.subgraph_search(K_n1, induced=True)
            sage: s_P = g.subgraph_search(P_n2, induced=True)
            sage: K_n1.is_isomorphic(s_K)
            True
            sage: P_n2.is_isomorphic(s_P)
            True

        Create several barbell graphs in a Sage graphics array::

            sage: g = []
            sage: j = []
            sage: for i in range(6):
            ...       k = graphs.BarbellGraph(i + 2, 4)
            ...       g.append(k)
            ...
            sage: for i in range(2):
            ...       n = []
            ...       for m in range(3):
            ...           n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...       j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time

        TESTS:

        The input ``n1`` must be `\geq 2`::

            sage: graphs.BarbellGraph(1, randint(0, 10^6))
            Traceback (most recent call last):
            ...
            ValueError: Invalid graph description, n1 should be >= 2
            sage: graphs.BarbellGraph(randint(-10^6, 1), randint(0, 10^6))
            Traceback (most recent call last):
            ...
            ValueError: Invalid graph description, n1 should be >= 2

        The input ``n2`` must be `\geq 0`::

            sage: graphs.BarbellGraph(randint(2, 10^6), -1)
            Traceback (most recent call last):
            ...
            ValueError: Invalid graph description, n2 should be >= 0
            sage: graphs.BarbellGraph(randint(2, 10^6), randint(-10^6, -1))
            Traceback (most recent call last):
            ...
            ValueError: Invalid graph description, n2 should be >= 0
            sage: graphs.BarbellGraph(randint(-10^6, 1), randint(-10^6, -1))
            Traceback (most recent call last):
            ...
            ValueError: Invalid graph description, n1 should be >= 2
        """
        # sanity checks
        if n1 < 2:
            raise ValueError("Invalid graph description, n1 should be >= 2")
        if n2 < 0:
            raise ValueError("Invalid graph description, n2 should be >= 0")

        pos_dict = {}

        for i in range(n1):
            x = float(cos((pi / 4) - ((2 * pi) / n1) * i) - (n2 / 2) - 1)
            y = float(sin((pi / 4) - ((2 * pi) / n1) * i) - (n2 / 2) - 1)
            j = n1 - 1 - i
            pos_dict[j] = (x, y)
        for i in range(n1, n1 + n2):
            x = float(i - n1 - (n2 / 2) + 1)
            y = float(i - n1 - (n2 / 2) + 1)
            pos_dict[i] = (x, y)
        for i in range(n1 + n2, (2 * n1) + n2):
            x = float(
                cos((5 * (pi / 4)) + ((2 * pi) / n1) * (i - n1 - n2))
                + (n2 / 2) + 2)
            y = float(
                sin((5 * (pi / 4)) + ((2 * pi) / n1) * (i - n1 - n2))
                + (n2 / 2) + 2)
            pos_dict[i] = (x, y)

        import networkx
        G = networkx.barbell_graph(n1, n2)
        return graph.Graph(G, pos=pos_dict, name="Barbell graph")


    def BuckyBall(self):
        r"""
        Create the Bucky Ball graph.

        This graph is a 3-regular 60-vertex planar graph. Its vertices
        and edges correspond precisely to the carbon atoms and bonds
        in buckminsterfullerene.  When embedded on a sphere, its 12
        pentagon and 20 hexagon faces are arranged exactly as the
        sections of a soccer ball.

        EXAMPLES:

        The Bucky Ball is planar. ::

            sage: g = graphs.BuckyBall()
            sage: g.is_planar()
            True

        The Bucky Ball can also be created by extracting the 1-skeleton
        of the Bucky Ball polyhedron, but this is much slower. ::

            sage: g = polytopes.buckyball().vertex_graph()
            sage: g.remove_loops()
            sage: h = graphs.BuckyBall()
            sage: g.is_isomorphic(h)
            True

        The graph is returned along with an attractive embedding. ::

            sage: g = graphs.BuckyBall()
            sage: g.plot(vertex_labels=False, vertex_size=10).show() # long time
        """
        edges = [(0, 2), (0, 48), (0, 59), (1, 3), (1, 9), (1, 58),
                 (2, 3), (2, 36), (3, 17), (4, 6), (4, 8), (4, 12),
                 (5, 7), (5, 9), (5, 16), (6, 7), (6, 20), (7, 21),
                 (8, 9), (8, 56), (10, 11), (10, 12), (10, 20), (11, 27),
                 (11, 47), (12, 13), (13, 46), (13, 54), (14, 15), (14, 16),
                 (14, 21), (15, 25), (15, 41), (16, 17), (17, 40), (18, 19),
                 (18, 20), (18, 26), (19, 21), (19, 24), (22, 23), (22, 31),
                 (22, 34), (23, 25), (23, 38), (24, 25), (24, 30), (26, 27),
                 (26, 30), (27, 29), (28, 29), (28, 31), (28, 35), (29, 44),
                 (30, 31), (32, 34), (32, 39), (32, 50), (33, 35), (33, 45),
                 (33, 51), (34, 35), (36, 37), (36, 40), (37, 39), (37, 52),
                 (38, 39), (38, 41), (40, 41), (42, 43), (42, 46), (42, 55),
                 (43, 45), (43, 53), (44, 45), (44, 47), (46, 47), (48, 49),
                 (48, 52), (49, 53), (49, 57), (50, 51), (50, 52), (51, 53),
                 (54, 55), (54, 56), (55, 57), (56, 58), (57, 59), (58, 59)
                 ]
        g = graph.Graph()
        g.add_edges(edges)
        g.name("Bucky Ball")

        pos = {
            0 :  (1.00000000000000, 0.000000000000000),
            1 :  (-1.00000000000000, 0.000000000000000),
            2 :  (0.500000000000000, 0.866025403784439),
            3 :  (-0.500000000000000, 0.866025403784439),
            4 :  (-0.252886764483159, -0.146004241548845),
            5 :  (-0.368953972399043, 0.0928336233191176),
            6 :  (-0.217853192651371, -0.0480798425451855),
            7 :  (-0.255589950938772, 0.0495517623332213),
            8 :  (-0.390242139418333, -0.225306404242310),
            9 :  (-0.586398703939125, -0.0441575936410641),
            10:  (-0.113926229169631, -0.101751920396670),
            11:  (-0.0461308635969359, -0.0928422349110366),
            12:  (-0.150564961379772, -0.164626477859040),
            13:  (-0.0848818904865275, -0.246123271631605),
            14:  (-0.170708060452244, 0.196571509298384),
            15:  (-0.0672882312715990, 0.212706320404226),
            16:  (-0.264873262319233, 0.273106701265196),
            17:  (-0.254957754106411, 0.529914971178085),
            18:  (-0.103469165775548, 0.00647061768205703),
            19:  (-0.113590051906687, 0.0655812470455896),
            20:  (-0.145082862532183, -0.0477870484199328),
            21:  (-0.179962687765901, 0.103901506225732),
            22:  (0.0573383021786124, 0.0863716172289798),
            23:  (0.0311566333625530, 0.149538968816603),
            24:  (-0.0573383021786121, 0.0863716172289799),
            25:  (-0.0311566333625527, 0.149538968816603),
            26:  (-0.0517345828877740, 0.00161765442051429),
            27:  (-0.0244663616211774, -0.0456122902452611),
            28:  (0.0517345828877743, 0.00161765442051431),
            29:  (0.0244663616211777, -0.0456122902452611),
            30:  (-0.0272682212665964, 0.0439946358247470),
            31:  (0.0272682212665968, 0.0439946358247470),
            32:  (0.179962687765901, 0.103901506225732),
            33:  (0.145082862532184, -0.0477870484199329),
            34:  (0.113590051906687, 0.0655812470455895),
            35:  (0.103469165775548, 0.00647061768205698),
            36:  (0.254957754106411, 0.529914971178085),
            37:  (0.264873262319233, 0.273106701265196),
            38:  (0.0672882312715993, 0.212706320404226),
            39:  (0.170708060452245, 0.196571509298384),
            40:  (1.59594559789866e-16, 0.450612808484620),
            41:  (2.01227923213310e-16, 0.292008483097691),
            42:  (0.0848818904865278, -0.246123271631605),
            43:  (0.150564961379773, -0.164626477859040),
            44:  (0.0461308635969362, -0.0928422349110366),
            45:  (0.113926229169631, -0.101751920396670),
            46:  (1.66533453693773e-16, -0.207803012451463),
            47:  (1.80411241501588e-16, -0.131162494091179),
            48:  (0.586398703939126, -0.0441575936410641),
            49:  (0.390242139418333, -0.225306404242310),
            50:  (0.255589950938772, 0.0495517623332212),
            51:  (0.217853192651372, -0.0480798425451855),
            52:  (0.368953972399044, 0.0928336233191175),
            53:  (0.252886764483159, -0.146004241548845),
            54:  (-0.104080710079810, -0.365940324584313),
            55:  (0.104080710079811, -0.365940324584313),
            56:  (-0.331440949832714, -0.485757377537020),
            57:  (0.331440949832715, -0.485757377537021),
            58:  (-0.500000000000000, -0.866025403784438),
            59:  (0.500000000000000, -0.866025403784439)
        }

        g.set_pos(pos)

        return g

    def BullGraph(self):
        r"""
        Returns a bull graph with 5 nodes.

        A bull graph is named for its shape. It's a triangle with horns.
        This constructor depends on `NetworkX <http://networkx.lanl.gov>`_
        numeric labeling. For more information, see this
        `Wikipedia article on the bull graph <http://en.wikipedia.org/wiki/Bull_graph>`_.

        PLOTTING:

        Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the bull graph
        is drawn as a triangle with the first node (0) on the bottom. The
        second and third nodes (1 and 2) complete the triangle. Node 3 is
        the horn connected to 1 and node 4 is the horn connected to node
        2.

        ALGORITHM:

        Uses `NetworkX <http://networkx.lanl.gov>`_.

        EXAMPLES:

        Construct and show a bull graph::

            sage: g = graphs.BullGraph(); g
            Bull graph: Graph on 5 vertices
            sage: g.show() # long time

        The bull graph has 5 vertices and 5 edges. Its radius is 2, its
        diameter 3, and its girth 3. The bull graph is planar with chromatic
        number 3 and chromatic index also 3. ::

            sage: g.order(); g.size()
            5
            5
            sage: g.radius(); g.diameter(); g.girth()
            2
            3
            3
            sage: g.chromatic_number()
            3

        The bull graph has chromatic polynomial `x(x - 2)(x - 1)^3` and
        Tutte polynomial `x^4 + x^3 + x^2 y`. Its characteristic polynomial
        is `x(x^2 - x - 3)(x^2 + x - 1)`, which follows from the definition of
        characteristic polynomials for graphs, i.e. `\det(xI - A)`, where
        `x` is a variable, `A` the adjacency matrix of the graph, and `I`
        the identity matrix of the same dimensions as `A`. ::

            sage: chrompoly = g.chromatic_polynomial()
            sage: bool(expand(x * (x - 2) * (x - 1)^3) == chrompoly)
            True
            sage: charpoly = g.characteristic_polynomial()
            sage: M = g.adjacency_matrix(); M
            [0 1 1 0 0]
            [1 0 1 1 0]
            [1 1 0 0 1]
            [0 1 0 0 0]
            [0 0 1 0 0]
            sage: Id = identity_matrix(ZZ, M.nrows())
            sage: D = x*Id - M
            sage: bool(D.determinant() == charpoly)
            True
            sage: bool(expand(x * (x^2 - x - 3) * (x^2 + x - 1)) == charpoly)
            True
        """
        pos_dict = {0:(0,0), 1:(-1,1), 2:(1,1), 3:(-2,2), 4:(2,2)}
        import networkx
        G = networkx.bull_graph()
        return graph.Graph(G, pos=pos_dict, name="Bull graph")

    def ButterflyGraph(self):
        r"""
        Returns the butterfly graph.

        Let `C_3` be the cycle graph on 3 vertices. The butterfly or bowtie
        graph is obtained by joining two copies of `C_3` at a common vertex,
        resulting in a graph that is isomorphic to the friendship graph `F_2`.
        For more information, see this
        `Wikipedia article on the butterfly graph <http://en.wikipedia.org/wiki/Butterfly_graph>`_.

        .. seealso::

            - :meth:`GraphGenerators.FriendshipGraph`

        EXAMPLES:

        The butterfly graph is a planar graph on 5 vertices and having
        6 edges. ::

            sage: G = graphs.ButterflyGraph(); G
            Butterfly graph: Graph on 5 vertices
            sage: G.show()  # long time
            sage: G.is_planar()
            True
            sage: G.order()
            5
            sage: G.size()
            6

        It has diameter 2, girth 3, and radius 1. ::

            sage: G.diameter()
            2
            sage: G.girth()
            3
            sage: G.radius()
            1

        The butterfly graph is Eulerian, with chromatic number 3. ::

            sage: G.is_eulerian()
            True
            sage: G.chromatic_number()
            3
        """
        edge_dict = {
            0: [3,4],
            1: [2,4],
            2: [4],
            3: [4]}
        pos_dict = {
            0: [-1, 1],
            1: [1, 1],
            2: [1, -1],
            3: [-1, -1],
            4: [0, 0]}
        return graph.Graph(edge_dict, pos=pos_dict, name="Butterfly graph")

    def CircularLadderGraph(self, n):
        """
        Returns a circular ladder graph with 2\*n nodes.

        A Circular ladder graph is a ladder graph that is connected at the
        ends, i.e.: a ladder bent around so that top meets bottom. Thus it
        can be described as two parallel cycle graphs connected at each
        corresponding node pair.

        This constructor depends on NetworkX numeric labels.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the circular
        ladder graph is displayed as an inner and outer cycle pair, with
        the first n nodes drawn on the inner circle. The first (0) node is
        drawn at the top of the inner-circle, moving clockwise after that.
        The outer circle is drawn with the (n+1)th node at the top, then
        counterclockwise as well.

        EXAMPLES: Construct and show a circular ladder graph with 26 nodes

        ::

            sage: g = graphs.CircularLadderGraph(13)
            sage: g.show() # long time

        Create several circular ladder graphs in a Sage graphics array

        ::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.CircularLadderGraph(i+3)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time
        """
        pos_dict = {}
        for i in range(n):
            x = float(cos((pi/2) + ((2*pi)/n)*i))
            y = float(sin((pi/2) + ((2*pi)/n)*i))
            pos_dict[i] = [x,y]
        for i in range(n,2*n):
            x = float(2*(cos((pi/2) + ((2*pi)/n)*(i-n))))
            y = float(2*(sin((pi/2) + ((2*pi)/n)*(i-n))))
            pos_dict[i] = (x,y)
        import networkx
        G = networkx.circular_ladder_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Circular Ladder graph")

    def ClawGraph(self):
        """
        Returns a claw graph.

        A claw graph is named for its shape. It is actually a complete
        bipartite graph with (n1, n2) = (1, 3).

        PLOTTING: See CompleteBipartiteGraph.

        EXAMPLES: Show a Claw graph

        ::

            sage: (graphs.ClawGraph()).show() # long time

        Inspect a Claw graph

        ::

            sage: G = graphs.ClawGraph()
            sage: G
            Claw graph: Graph on 4 vertices
        """
        pos_dict = {0:(0,1),1:(-1,0),2:(0,0),3:(1,0)}
        import networkx
        G = networkx.complete_bipartite_graph(1,3)
        return graph.Graph(G, pos=pos_dict, name="Claw graph")

    def CycleGraph(self, n):
        r"""
        Returns a cycle graph with n nodes.

        A cycle graph is a basic structure which is also typically called
        an n-gon.

        This constructor is dependent on vertices numbered 0 through n-1 in
        NetworkX ``cycle_graph()``

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, each cycle
        graph will be displayed with the first (0) node at the top, with
        the rest following in a counterclockwise manner.

        The cycle graph is a good opportunity to compare efficiency of
        filling a position dictionary vs. using the spring-layout algorithm
        for plotting. Because the cycle graph is very symmetric, the
        resulting plots should be similar (in cases of small n).

        Filling the position dictionary in advance adds O(n) to the
        constructor.

        EXAMPLES: Compare plotting using the predefined layout and
        networkx::

            sage: import networkx
            sage: n = networkx.cycle_graph(23)
            sage: spring23 = Graph(n)
            sage: posdict23 = graphs.CycleGraph(23)
            sage: spring23.show() # long time
            sage: posdict23.show() # long time

        We next view many cycle graphs as a Sage graphics array. First we
        use the ``CycleGraph`` constructor, which fills in the
        position dictionary::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.CycleGraph(i+3)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time

        Compare to plotting with the spring-layout algorithm::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    spr = networkx.cycle_graph(i+3)
            ...    k = Graph(spr)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time
        """
        pos_dict = {}
        for i in range(n):
            x = float(cos((pi/2) + ((2*pi)/n)*i))
            y = float(sin((pi/2) + ((2*pi)/n)*i))
            pos_dict[i] = (x,y)
        import networkx
        G = networkx.cycle_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Cycle graph")

    def DiamondGraph(self):
        """
        Returns a diamond graph with 4 nodes.

        A diamond graph is a square with one pair of diagonal nodes
        connected.

        This constructor depends on NetworkX numeric labeling.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the diamond
        graph is drawn as a diamond, with the first node on top, second on
        the left, third on the right, and fourth on the bottom; with the
        second and third node connected.

        EXAMPLES: Construct and show a diamond graph

        ::

            sage: g = graphs.DiamondGraph()
            sage: g.show() # long time
        """
        pos_dict = {0:(0,1),1:(-1,0),2:(1,0),3:(0,-1)}
        import networkx
        G = networkx.diamond_graph()
        return graph.Graph(G, pos=pos_dict, name="Diamond Graph")

    def EmptyGraph(self):
        """
        Returns an empty graph (0 nodes and 0 edges).

        This is useful for constructing graphs by adding edges and vertices
        individually or in a loop.

        PLOTTING: When plotting, this graph will use the default
        spring-layout algorithm, unless a position dictionary is
        specified.

        EXAMPLES: Add one vertex to an empty graph and then show::

            sage: empty1 = graphs.EmptyGraph()
            sage: empty1.add_vertex()
            0
            sage: empty1.show() # long time

        Use for loops to build a graph from an empty graph::

            sage: empty2 = graphs.EmptyGraph()
            sage: for i in range(5):
            ...    empty2.add_vertex() # add 5 nodes, labeled 0-4
            ...
            0
            1
            2
            3
            4
            sage: for i in range(3):
            ...    empty2.add_edge(i,i+1) # add edges {[0:1],[1:2],[2:3]}
            ...
            sage: for i in range(4)[1:]:
            ...    empty2.add_edge(4,i) # add edges {[1:4],[2:4],[3:4]}
            ...
            sage: empty2.show() # long time
        """
        return graph.Graph(sparse=True)

    def ToroidalGrid2dGraph(self,n1,n2):
        r"""
        Returns a toroidal 2-dimensional grid graph with `n_1n_2` nodes
        (`n_1` rows and `n_2` columns).

        The toroidal 2-dimensional grid with parameters `n_1,n_2` is
        the 2-dimensional grid graph with identical parameters
        to which are added the edges `((i,0),(i,n_2-1))` and
        `((0,i),(n_1-1,i))`.

        EXAMPLE:

        The toroidal 2-dimensional grid is a regular graph, while
        the usual 2-dimensional grid is not ::

            sage: tgrid = graphs.ToroidalGrid2dGraph(8,9)
            sage: print tgrid
            Toroidal 2D Grid Graph with parameters 8,9
            sage: grid = graphs.Grid2dGraph(8,9)
            sage: grid.is_regular()
            False
            sage: tgrid.is_regular()
            True
        """

        g = self.Grid2dGraph(n1,n2)

        g.add_edges([((i,0),(i,n2-1)) for i in range(n1)] + [((0,i),(n1-1,i)) for i in range(n2)])

        g.name("Toroidal 2D Grid Graph with parameters "+str(n1)+","+str(n2))

        return g

    def Grid2dGraph(self, n1, n2):
        r"""
        Returns a `2`-dimensional grid graph with `n_1n_2` nodes (`n_1` rows and
        `n_2` columns).

        A 2d grid graph resembles a `2` dimensional grid. All inner nodes are
        connected to their `4` neighbors. Outer (non-corner) nodes are
        connected to their `3` neighbors. Corner nodes are connected to their
        2 neighbors.

        This constructor depends on NetworkX numeric labels.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, nodes are
        labelled in (row, column) pairs with `(0, 0)` in the top left corner.
        Edges will always be horizontal and vertical - another advantage of
        filling the position dictionary.

        EXAMPLES: Construct and show a grid 2d graph Rows = `5`, Columns = `7`

        ::

            sage: g = graphs.Grid2dGraph(5,7)
            sage: g.show() # long time
        """
        pos_dict = {}
        for i in range(n1):
            y = -i
            for j in range(n2):
                x = j
                pos_dict[i,j] = (x,y)
        import networkx
        G = networkx.grid_2d_graph(n1,n2)
        return graph.Graph(G, pos=pos_dict, name="2D Grid Graph")

    def GridGraph(self, dim_list):
        """
        Returns an n-dimensional grid graph.

        INPUT:


        -  ``dim_list`` - a list of integers representing the
           number of nodes to extend in each dimension.


        PLOTTING: When plotting, this graph will use the default
        spring-layout algorithm, unless a position dictionary is
        specified.

        EXAMPLES::

            sage: G = graphs.GridGraph([2,3,4])
            sage: G.show()  # long time

        ::

            sage: C = graphs.CubeGraph(4)
            sage: G = graphs.GridGraph([2,2,2,2])
            sage: C.show()  # long time
            sage: G.show()  # long time
        """
        import networkx
        dim = [int(a) for a in dim_list]
        G = networkx.grid_graph(dim)
        return graph.Graph(G, name="Grid Graph for %s"%dim)

    def HanoiTowerGraph(self, pegs, disks, labels=True, positions=True):
        r"""
        Returns the graph whose vertices are the states of the
        Tower of Hanoi puzzle, with edges representing legal moves between states.

        INPUT:

        - ``pegs`` - the number of pegs in the puzzle, 2 or greater
        - ``disks`` - the number of disks in the puzzle, 1 or greater
        - ``labels`` - default: ``True``, if ``True`` the graph contains
          more meaningful labels, see explanation below.  For large instances,
          turn off labels for much faster creation of the graph.
        - ``positions`` - default: ``True``, if ``True`` the graph contains
          layout information.  This creates a planar layout for the case
          of three pegs.  For large instances, turn off layout information
          for much faster creation of the graph.

        OUTPUT:

        The Tower of Hanoi puzzle has a certain number of identical pegs
        and a certain number of disks, each of a different radius.
        Initially the disks are all on a single peg, arranged
        in order of their radii, with the largest on the bottom.

        The goal of the puzzle is to move the disks to any other peg,
        arranged in the same order.  The one constraint is that the
        disks resident on any one peg must always be arranged with larger
        radii lower down.

        The vertices of this graph represent all the possible states
        of this puzzle.  Each state of the puzzle is a tuple with length
        equal to the number of disks, ordered by largest disk first.
        The entry of the tuple is the peg where that disk resides.
        Since disks on a given peg must go down in size as we go
        up the peg, this totally describes the state of the puzzle.

        For example ``(2,0,0)`` means the large disk is on peg 2, the
        medium disk is on peg 0, and the small disk is on peg 0
        (and we know the small disk must be above the medium disk).
        We encode these tuples as integers with a base equal to
        the number of pegs, and low-order digits to the right.

        Two vertices are adjacent if we can change the puzzle from
        one state to the other by moving a single disk.  For example,
        ``(2,0,0)`` is adjacent to ``(2,0,1)`` since we can move
        the small disk off peg 0 and onto (the empty) peg 1.
        So the solution to a 3-disk puzzle (with at least
        two pegs) can be expressed by the shortest path between
        ``(0,0,0)`` and ``(1,1,1)``.  For more on this representation
        of the graph, or its properties, see [ARETT-DOREE]_.

        For greatest speed we create graphs with integer vertices,
        where we encode the tuples as integers with a base equal
        to the number of pegs, and low-order digits to the right.
        So for example, in a 3-peg puzzle with 5 disks, the
        state ``(1,2,0,1,1)`` is encoded as
        `1\ast 3^4 + 2\ast 3^3 + 0\ast 3^2 + 1\ast 3^1 + 1\ast 3^0 = 139`.

        For smaller graphs, the labels that are the tuples are informative,
        but slow down creation of the graph.  Likewise computing layout
        information also incurs a significant speed penalty. For maximum
        speed, turn off labels and layout and decode the
        vertices explicitly as needed.  The
        :meth:`sage.rings.integer.Integer.digits`
        with the ``padsto`` option is a quick way to do this, though you
        may want to reverse the list that is output.

        PLOTTING:

        The layout computed when ``positions = True`` will
        look especially good for the three-peg case, when the graph is known
        to be planar.  Except for two small cases on 4 pegs, the graph is
        otherwise not planar, and likely there is a better way to layout
        the vertices.

        EXAMPLES:

        A classic puzzle uses 3 pegs.  We solve the 5 disk puzzle using
        integer labels and report the minimum number of moves required.
        Note that `3^5-1` is the state where all 5 disks
        are on peg 2. ::

            sage: H = graphs.HanoiTowerGraph(3, 5, labels=False, positions=False)
            sage: H.distance(0, 3^5-1)
            31

        A slightly larger instance. ::

            sage: H = graphs.HanoiTowerGraph(4, 6, labels=False, positions=False)
            sage: H.num_verts()
            4096
            sage: H.distance(0, 4^6-1)
            17

        For a small graph, labels and layout information can be useful.
        Here we explicitly list a solution as a list of states. ::

            sage: H = graphs.HanoiTowerGraph(3, 3, labels=True, positions=True)
            sage: H.shortest_path((0,0,0), (1,1,1))
            [(0, 0, 0), (0, 0, 1), (0, 2, 1), (0, 2, 2), (1, 2, 2), (1, 2, 0), (1, 1, 0), (1, 1, 1)]

        Some facts about this graph with `p` pegs and `d` disks:

        - only automorphisms are the "obvious" ones - renumber the pegs.
        - chromatic number is less than or equal to `p`
        - independence number is `p^{d-1}`

        ::

            sage: H = graphs.HanoiTowerGraph(3,4,labels=False,positions=False)
            sage: H.automorphism_group().is_isomorphic(SymmetricGroup(3))
            True
            sage: H.chromatic_number()
            3
            sage: len(H.independent_set()) == 3^(4-1)
            True

        TESTS:

        It is an error to have just one peg (or less). ::

            sage: graphs.HanoiTowerGraph(1, 5)
            Traceback (most recent call last):
            ...
            ValueError: Pegs for Tower of Hanoi graph should be two or greater (not 1)

        It is an error to have zero disks (or less). ::

            sage: graphs.HanoiTowerGraph(2, 0)
            Traceback (most recent call last):
            ...
            ValueError: Disks for Tower of Hanoi graph should be one or greater (not 0)

        .. rubric:: Citations

        .. [ARETT-DOREE] Arett, Danielle and Doree, Suzanne
           "Coloring and counting on the Hanoi graphs"
           Mathematics Magazine, Volume 83, Number 3, June 2010, pages 200-9


        AUTHOR:

        - Rob Beezer, (2009-12-26), with assistance from Su Doree

        """

        # sanitize input
        from sage.rings.all import Integer
        pegs = Integer(pegs)
        if pegs < 2:
            raise ValueError("Pegs for Tower of Hanoi graph should be two or greater (not %d)" % pegs)
        disks = Integer(disks)
        if disks < 1:
            raise ValueError("Disks for Tower of Hanoi graph should be one or greater (not %d)" % disks)

        # Each state of the puzzle is a tuple with length
        # equal to the number of disks, ordered by largest disk first
        # The entry of the tuple is the peg where that disk resides
        # Since disks on a given peg must go down in size as we go
        # up the peg, this totally describes the puzzle
        # We encode these tuples as integers with a base equal to
        # the number of pegs, and low-order digits to the right

        # complete graph on number of pegs when just a single disk
        edges = [[i,j] for i in range(pegs) for j in range(i+1,pegs)]

        nverts = 1
        for d in range(2, disks+1):
            prevedges = edges      # remember subgraph to build from
            nverts = pegs*nverts   # pegs^(d-1)
            edges = []

            # Take an edge, change its two states in the same way by adding
            # a large disk to the bottom of the same peg in each state
            # This is accomplished by adding a multiple of pegs^(d-1)
            for p in range(pegs):
                largedisk = p*nverts
                for anedge in prevedges:
                    edges.append([anedge[0]+largedisk, anedge[1]+largedisk])

            # Two new states may only differ in the large disk
            # being the only disk on two different pegs, thus
            # otherwise being a common state with one less disk
            # We construct all such pairs of new states and add as edges
            from sage.combinat.subset import Subsets
            for state in range(nverts):
                emptypegs = range(pegs)
                reduced_state = state
                for i in range(d-1):
                    apeg = reduced_state % pegs
                    if apeg in emptypegs:
                        emptypegs.remove(apeg)
                    reduced_state = reduced_state//pegs
                for freea, freeb in Subsets(emptypegs, 2):
                    edges.append([freea*nverts+state,freeb*nverts+state])

        H = graph.Graph({}, loops=False, multiedges=False)
        H.add_edges(edges)


        # Making labels and/or computing positions can take a long time,
        # relative to just constructing the edges on integer vertices.
        # We try to minimize coercion overhead, but need Sage
        # Integers in order to use digits() for labels.
        # Getting the digits with custom code was no faster.
        # Layouts are circular (symmetric on the number of pegs)
        # radiating outward to the number of disks (radius)
        # Algorithm uses some combination of alternate
        # clockwise/counterclockwise placements, which
        # works well for three pegs (planar layout)
        #
        from sage.functions.trig import sin, cos, csc
        if labels or positions:
            mapping = {}
            pos = {}
            a = Integer(-1)
            one = Integer(1)
            if positions:
                radius_multiplier = 1 + csc(pi/pegs)
                sine = []; cosine = []
                for i in range(pegs):
                    angle = 2*i*pi/float(pegs)
                    sine.append(sin(angle))
                    cosine.append(cos(angle))
            for i in range(pegs**disks):
                a += one
                state = a.digits(base=pegs, padto=disks)
                if labels:
                    state.reverse()
                    mapping[i] = tuple(state)
                    state.reverse()
                if positions:
                    locx = 0.0; locy = 0.0
                    radius = 1.0
                    parity = -1.0
                    for index in range(disks):
                        p = state[index]
                        radius *= radius_multiplier
                        parity *= -1.0
                        locx_temp = cosine[p]*locx - parity*sine[p]*locy + radius*cosine[p]
                        locy_temp = parity*sine[p]*locx + cosine[p]*locy - radius*parity*sine[p]
                        locx = locx_temp
                        locy = locy_temp
                    pos[i] = (locx,locy)
            # set positions, then relabel (not vice versa)
            if positions:
                H.set_pos(pos)
            if labels:
                H.relabel(mapping)

        return H

    def HallJankoGraph(self, from_string=True):
        r"""
        Returns the Hall-Janko graph.

        For more information on the Hall-Janko graph, see its
        :wikipedia:`Wikipedia page <Hall-Janko_graph>`.

        The construction used to generate this graph in Sage is by
        a 100-point permutation representation of the Janko group `J_2`,
        as described in version 3 of the ATLAS of Finite Group
        representations, in particular on the page `ATLAS: J2
        — Permutation representation on 100 points
        <http://brauer.maths.qmul.ac.uk/Atlas/v3/permrep/J2G1-p100B0>`_.

        INPUT:

        - ``from_string`` (boolean) -- whether to build the graph from
          its sparse6 string or through GAP. The two methods return the
          same graph though doing it through GAP takes more time. It is
          set to ``True`` by default.

        EXAMPLES::

            sage: g = graphs.HallJankoGraph()
            sage: g.is_regular(36)
            True
            sage: g.is_vertex_transitive()
            True

        Is it really strongly regular with parameters 14, 12? ::

            sage: nu = set(g.neighbors(0))
            sage: for v in range(1, 100):
            ...     if v in nu:
            ...        expected = 14
            ...     else:
            ...        expected = 12
            ...     nv = set(g.neighbors(v))
            ...     nv.discard(0)
            ...     if len(nu & nv) != expected:
            ...        print "Something is wrong here!!!"
            ...        break

        Some other properties that we know how to check::

            sage: g.diameter()
            2
            sage: g.girth()
            3
            sage: factor(g.characteristic_polynomial())
            (x - 36) * (x - 6)^36 * (x + 4)^63

        TESTS::

            sage: gg = graphs.HallJankoGraph(from_string=False) # long time
            sage: g == gg # long time
            True
        """

        string = (":~?@c__E@?g?A?w?A@GCA_?CA`OWF`W?EAW?@?_OD@_[GAgcIaGGB@OcIA"
                  "wCE@o_K_?GB@?WGAouC@OsN_?GB@O[GB`A@@_e?@OgLB_{Q_?GC@O[GAOs"
                  "OCWGBA?kKBPA@?_[KB_{OCPKT`o_RD`]A?o[HBOwODW?DA?cIB?wRDP[X`"
                  "ogKB_{QD@]B@o_KBPWXE`mC@o_JB?{PDPq@?oWGA_{OCPKTDp_YEwCA@_c"
                  "IBOwOC`OX_OGB@?WPDPcYFg?C@_gKBp?SE@cYF`{_`?SGAOoOC`_\\FwCE"
                  "A?gKBO{QD@k[FqI??_OFA_oQE@k\\Fq?`GgCB@pGRD@_XFP{a_?SE@ocIA"
                  "ooNCPOUEqU@?oODA?cJB_{UEqYC@_kLC@CREPk]GAGbHgCA@?SMBpCSD`["
                  "YFq?`Ga]BA?gPC`KSD`_\\Fa?cHWGB@?[IAooPD`[WF@s^HASeIg?@@OcP"
                  "C`KYF@w^GQ[h`O[HAooMC@CQCpSVEPk\\GaSeIG?FA?kLB_{OC`OVE@cYG"
                  "QUA@?WLBp?PC`KVEqKgJg?DA?sMBpCSDP[WEQKfIay@?_KD@_[GC`SUE@k"
                  "[FaKdHa[k_?OLC@CRD@WVEpo^HAWfIAciIqoo_?CB@?kMCpOUE`o\\GAKg"
                  "IQgq_?GD@_[GB?{OCpWVE@cYFACaHAWhJR?q_?CC@_kKBpC\\GACdHa[kJ"
                  "a{o_?CA?oOFBpGRD@o\\GaKdIQonKrOt_?WHA`?PC`KTD`k]FqSeIaolJr"
                  "CqLWCA@OkKCPGRDpcYGAKdIAgjJAsmJr?t__OE@ogJB_{XEps`HA[gIQwn"
                  "KWKGAOoMBpGUE`k[Fa?aHqckJbSuLw?@?_SHA_kLC@OTFPw^GaOkLg?B@?"
                  "[HA_{PDP_XFaCbHa[gIqooKRWx_?CFBpOTE@cZFPw^GACcHQgoKrSvMwWG"
                  "BOwQCp_YFP{`HASfJAwnKRSx_OSSDP[WEq?aGqSfIQsoKR_zNWCE@o_HA_"
                  "sREPg^GAGcHQWfIAciKbOxNg?A@__IAooMC`KTD`g\\GAKcIasoKrOtLb["
                  "wMbyCA?cKBp?TD`[WE`s^GQGbHqcjJrK{NRw~_oODA?sNC@CQCpOZF@s]G"
                  "QOfIaolJrGsLbk}_?OFA_sRD@SVE`k[HQcjJa{qLb[xMb|?_OOFA?cIAos"
                  "RDP_ZFa?aGqOfIAsuMbk{Ns@@OsQAA_sPDPWXE`o\\FqKdIQkkJrCuLr_x"
                  "Mro}NsDAPG?@@OWFApKUE@o`IQolKRKsLrc|NsQC@OWGAOgJCpOWE`o_GQ"
                  "KiIqwnKr_~OcLCPS]A?oWHA_oMBpKSDP[\\FagjKBWxMbk{OSQ@@O_IAoo"
                  "LBpCSD`g\\FaGbHQWgIQgmKRKwMRl?PgGC@OWHB@KSE@c[FqCaGqSeIAkk"
                  "KBCqLBSuMBpGQWCA@?cKBOwRDPWVE@k^GqOfJr?pKbKtLrs}OSHDQwKIBO"
                  "wPD@WWEQ?`HQWfIQglKBOtLbo}Ns@@OsTE_?kLCpWWHA[gIqomKBGwMRgz"
                  "NBw~OSPDPc\\H_?CFAOoLCPSVE`o\\GAOeJAwpKbKtMrx?Qcq??OKFA?gJ"
                  "B`?QDpcYEpo]FqKfIAgjJB?qKr_{NS@A__SE@o_HBO{PC`OTD`{_HaciIq"
                  "{vMbt?OcPFQCeB@?SKBOwRD@SXE`k[FPw`HQ_lKRKxNRxBPC\\HQclK_?K"
                  "EB?sOC`OTDa?`GqWgJRCrNBw~OSHFQStMRtDQ_?KC@OoQE`k_GaOdHa[gI"
                  "q{tMBg|Nb|?OcPMSDDQSwCB@_cJB_{OCpOVFP{dHa[jJQwqKrk}NsHBQCd"
                  "MRtMA?oSEA_wPDp_YEpo]GAOeIq{pLBk}NsLEQCtNTDU??OKEA_oLC@[[G"
                  "aKnKBOtLbk~OCPFQStNSDLSTgGKC@GSD`[WEpw_GQGcIAciJAwpKb_xMbk"
                  "~QShJRc|R`_wNCPcZF@s^GAGbHA_hJR?qKrOvMRg|NsDEPsxTTgCB@?gJB"
                  "?sMC@CUDp_]FqCaHQcjJQwtLrhCPS\\IRCtQTw?B@?SHA_wPC`_aGqOiJa"
                  "{oKRKvMRpFQChKRtXVUTi??ocNC@KUE@cYFaGdHa_mJrKsLb[yMro|OcXI"
                  "RdPTTddZaOgJB@?UEPk[FQCfIaolJrSvMBczNR|AOsXFQCtOTtaB@?WGAP"
                  "?TEPo\\GAGdHqgmKBCqLR[xMb|?PC`HQs|TTt`XUtu@?o[HB?sNCPGXF@{"
                  "_GQKcIqolJb_yNCLDPs`MRtDRTTdYUwSEA?kLB`CWF@s]FqGgIqooLRgzN"
                  "RxFQSlMSDDQTDXVUTi@?_KDAOoLBpKUEQOfIa{oLB_xMrt?Os\\HQcpMST"
                  "HSTtl[VT}A@ocJBOwSD`_XEpo_Ha_mJrKtLbgzNSTGQspLRtDUUDp\\WG["
                  "HB`CQCp[WFQGgIQgkJQ{rLbc{Nc@APsdLRt@PSt\\WUtt_Wn")

        if from_string:
            g = graph.Graph(string, loops = False, multiedges = False)
        else:

            # The following construction is due to version 3 of the ATLAS of
            # Finite Group Representations, specifically the page at
            # http://brauer.maths.qmul.ac.uk/Atlas/v3/permrep/J2G1-p100B0 .

            from sage.interfaces.gap import gap
            gap.eval("g1 := (1,84)(2,20)(3,48)(4,56)(5,82)(6,67)(7,55)(8,41)"
                     "(9,35)(10,40)(11,78)(12,100)(13,49)(14,37)(15,94)(16,76)"
                     "(17,19)(18,44)(21,34)(22,85)(23,92)(24,57)(25,75)(26,28)"
                     "(27,64)(29,90)(30,97)(31,38)(32,68)(33,69)(36,53)(39,61)"
                     "(42,73)(43,91)(45,86)(46,81)(47,89)(50,93)(51,96)(52,72)"
                     "(54,74)(58,99)(59,95)(60,63)(62,83)(65,70)(66,88)(71,87)"
                     "(77,98)(79,80);")

            gap.eval("g2 := (1,80,22)(2,9,11)(3,53,87)(4,23,78)(5,51,18)"
                     "(6,37,24)(8,27,60)(10,62,47)(12,65,31)(13,64,19)"
                     "(14,61,52)(15,98,25)(16,73,32)(17,39,33)(20,97,58)"
                     "(21,96,67)(26,93,99)(28,57,35)(29,71,55)(30,69,45)"
                     "(34,86,82)(38,59,94)(40,43,91)(42,68,44)(46,85,89)"
                     "(48,76,90)(49,92,77)(50,66,88)(54,95,56)(63,74,72)"
                     "(70,81,75)(79,100,83);")

            gap.eval("G := Group([g1,g2]);")
            edges = gap('Orbit(G,[1,5],OnSets)').sage()
            g = graph.Graph([(int(u), int(v)) for u,v in edges])
            g.relabel()

        _circle_embedding(g, range(100))
        g.name("Hall-Janko graph")
        return g

    def HararyGraph( self, k, n ):
        r"""
        Returns the Harary graph on `n` vertices and connectivity `k`, where
        `2 \leq k < n`.

        A `k`-connected graph `G` on `n` vertices requires the minimum degree
        `\delta(G)\geq k`, so the minimum number of edges `G` should have is
        `\lceil kn/2\rceil`. Harary graphs achieve this lower bound, that is,
        Harary graphs are minimal `k`-connected graphs on `n` vertices.

        The construction provided uses the method CirculantGraph.  For more
        details, see the book D. B. West, Introduction to Graph Theory, 2nd
        Edition, Prentice Hall, 2001, p. 150--151; or the `MathWorld article on
        Harary graphs <http://mathworld.wolfram.com/HararyGraph.html>`_.

        EXAMPLES:

        Harary graphs `H_{k,n}`::

            sage: h = graphs.HararyGraph(5,9); h
            Harary graph 5, 9: Graph on 9 vertices
            sage: h.order()
            9
            sage: h.size()
            23
            sage: h.vertex_connectivity()
            5

        TESTS:

        Connectivity of some Harary graphs::

            sage: n=10
            sage: for k in range(2,n):
            ...       g = graphs.HararyGraph(k,n)
            ...       if k != g.vertex_connectivity():
            ...          print "Connectivity of Harary graphs not satisfied."
        """
        if k < 2:
            raise ValueError("Connectivity parameter k should be at least 2.")
        if k >= n:
            raise ValueError("Number of vertices n should be greater than k.")

        if k%2 == 0:
            G = self.CirculantGraph( n, range(1,k/2+1) )
        else:
            if n%2 == 0:
                G = self.CirculantGraph( n, range(1,(k-1)/2+1) )
                for i in range(n):
                    G.add_edge( i, (i+n/2)%n )
            else:
                G = self.HararyGraph( k-1, n )
                for i in range((n-1)/2+1):
                    G.add_edge( i, (i+(n-1)/2)%n )
        G.name('Harary graph {0}, {1}'.format(k,n))
        return G

    def HarriesGraph(self, embedding=1):
        r"""
        Returns the Harries Graph.

        The Harries graph is a Hamiltonian 3-regular graph on 70
        vertices. See the :wikipedia:`Wikipedia page on the Harries
        graph <Harries_graph>`.

        The default embedding here is to emphasize the graph's 4 orbits.
        This graph actually has a funny construction. The following
        procedure gives an idea of it, though not all the adjacencies
        are being properly defined.

        #. Take two disjoint copies of a :meth:`Petersen graph
           <PetersenGraph>`. Their vertices will form an orbit of the
           final graph.

        #. Subdivide all the edges once, to create 15+15=30 new
           vertices, which together form another orbit.

        #. Create 15 vertices, each of them linked to 2 corresponding
           vertices of the previous orbit, one in each of the two
           subdivided Petersen graphs. At the end of this step all
           vertices from the previous orbit have degree 3, and the only
           vertices of degree 2 in the graph are those that were just
           created.

        #. Create 5 vertices connected only to the ones from the
           previous orbit so that the graph becomes 3-regular.

        INPUT:

        - ``embedding`` -- two embeddings are available, and can be
          selected by setting ``embedding`` to 1 or 2.

        EXAMPLES::

            sage: g = graphs.HarriesGraph()
            sage: g.order()
            70
            sage: g.size()
            105
            sage: g.girth()
            10
            sage: g.diameter()
            6
            sage: g.show(figsize=[10, 10])   # long time
            sage: graphs.HarriesGraph(embedding=2).show(figsize=[10, 10])   # long time

        TESTS::

            sage: graphs.HarriesGraph(embedding=3)
            Traceback (most recent call last):
            ...
            ValueError: The value of embedding must be 1 or 2.

        """
        g = graphs.LCFGraph(70, [-29, -19, -13, 13, 21, -27, 27, 33, -13, 13,
                                 19, -21, -33, 29], 5)
        g.name("Harries Graph")

        if embedding == 1:
            gpos = g.get_pos()
            ppos = graphs.PetersenGraph().get_pos()

            # The graph's four orbits
            o = [None]*4
            o[0] = [0, 2, 6, 8, 14, 16, 20, 22, 28, 30, 34, 36, 42, 44, 48, 50,
                    56, 58, 62, 64]
            o[1] = [1, 3, 5, 7, 9, 13, 15, 17, 19, 21, 23, 27, 29, 31, 33, 35,
                    37, 41, 43, 45, 47, 49, 51, 55, 57, 59, 61, 63, 65, 69]
            o[2] = [60, 10, 12, 4, 24, 26, 18, 38, 40, 32, 52, 54, 46, 66, 68]
            o[3] = [11, 25, 39, 53, 67]

            # Correspondence between the vertices of one of the two Petersen
            # graphs on o[0] and the vertices of a standard Petersen graph
            # object
            g_to_p = {0: 0, 2: 1, 42: 5, 44: 8, 14: 7, 16: 2, 56: 9, 58: 6,
                      28: 4, 30: 3}

            # Correspondence between the vertices of the other Petersen graph
            # on o[0] and the vertices of the first one
            g_to_g = {64: 44, 34: 0, 36: 28, 6: 2, 8: 58, 48: 16, 50: 30,
                      20: 14, 22: 56, 62: 42}

            # Position for the vertices from the first copy
            for v, i in g_to_p.iteritems():
                gpos[v] = ppos[i]

            # Position for the vertices in the second copy. Moves the first,
            # too.
            offset = 3.5
            for v, i in g_to_g.iteritems():
                x, y = gpos[i]
                gpos[v] = (x + offset*0.5, y)
                gpos[i] = (x - offset*0.5, y)

            # Vertices from o[1]. These are actually the "edges" of the
            # copies of Petersen.
            for v in o[1]:
                p1, p2 = [gpos[x] for x in g.neighbors(v) if x in o[0]]
                gpos[v] = ((p1[0] + p2[0])/2, (p1[1] + p2[1])/2)

            # 15 vertices from o[2]
            for i, v in enumerate(o[2]):
                gpos[v] = (-1.75 + i*.25, 2)

            # 5 vertices from o[3]
            for i, v in enumerate(o[3]):
                gpos[v] = (-1 + i*.5, 2.5)

            return g

        elif embedding == 2:
            return g
        else:
            raise ValueError("The value of embedding must be 1 or 2.")

    def HarriesWongGraph(self, embedding=1):
        r"""
        Returns the Harries-Wong Graph.

        See the :wikipedia:`Wikipedia page on the Harries-Wong graph
        <Harries-Wong_graph>`.

        *About the default embedding:*

        The default embedding is an attempt to emphasize the graph's
        8 (!!!) different orbits. In order to understand this better,
        one can picture the graph as being built in the following way:

            #. One first creates a 3-dimensional cube (8 vertices, 12
               edges), whose vertices define the first orbit of the
               final graph.

            #. The edges of this graph are subdivided once, to create 12
               new vertices which define a second orbit.

            #. The edges of the graph are subdivided once more, to
               create 24 new vertices giving a third orbit.

            #. 4 vertices are created and made adjacent to the vertices
               of the second orbit so that they have degree
               3. These 4 vertices also define a new orbit.

            #. In order to make the vertices from the third orbit
               3-regular (they all miss one edge), one creates a binary
               tree on 1 + 3 + 6 + 12 vertices. The leaves of this new
               tree are made adjacent to the 12 vertices of the third
               orbit, and the graph is now 3-regular. This binary tree
               contributes 4 new orbits to the Harries-Wong graph.

        INPUT:

        - ``embedding`` -- two embeddings are available, and can be
          selected by setting ``embedding`` to 1 or 2.

        EXAMPLES::

            sage: g = graphs.HarriesWongGraph()
            sage: g.order()
            70
            sage: g.size()
            105
            sage: g.girth()
            10
            sage: g.diameter()
            6
            sage: orbits = g.automorphism_group(orbits=True)[-1]
            sage: g.show(figsize=[15, 15], partition=orbits)   # long time

        Alternative embedding::

            sage: graphs.HarriesWongGraph(embedding=2).show()

        TESTS::

            sage: graphs.HarriesWongGraph(embedding=3)
            Traceback (most recent call last):
            ...
            ValueError: The value of embedding must be 1 or 2.
        """

        L = [9, 25, 31, -17, 17, 33, 9, -29, -15, -9, 9, 25, -25, 29, 17, -9,
             9, -27, 35, -9, 9, -17, 21, 27, -29, -9, -25, 13, 19, -9, -33,
             -17, 19, -31, 27, 11, -25, 29, -33, 13, -13, 21, -29, -21, 25,
             9, -11, -19, 29, 9, -27, -19, -13, -35, -9, 9, 17, 25, -9, 9, 27,
             -27, -21, 15, -9, 29, -29, 33, -9, -25]

        g = graphs.LCFGraph(70, L, 1)
        g.name("Harries-Wong graph")

        if embedding == 1:
            d = g.get_pos()

            # Binary tree (left side)
            d[66] = (-9.5, 0)
            _line_embedding(g, [37, 65, 67], first=(-8, 2.25),
                    last=(-8, -2.25))
            _line_embedding(g, [36, 38, 64, 24, 68, 30], first=(-7, 3),
                    last=(-7, -3))
            _line_embedding(g, [35, 39, 63, 25, 59, 29, 11, 5, 55, 23, 69, 31],
                    first=(-6, 3.5), last=(-6, -3.5))

            # Cube, corners: [9, 15, 21, 27, 45, 51, 57, 61]
            _circle_embedding(g, [61, 9], center=(0, -1.5), shift=.2,
                    radius=4)
            _circle_embedding(g, [27, 15], center=(0, -1.5), shift=.7,
                    radius=4*.707)
            _circle_embedding(g, [51, 21], center=(0, 2.5), shift=.2,
                    radius=4)
            _circle_embedding(g, [45, 57], center=(0, 2.5), shift=.7,
                    radius=4*.707)

            # Cube, subdivision
            _line_embedding(g, [21, 22, 43, 44, 45], first=d[21], last=d[45])
            _line_embedding(g, [21, 4, 3, 56, 57], first=d[21], last=d[57])
            _line_embedding(g, [57, 12, 13, 14, 15], first=d[57], last=d[15])
            _line_embedding(g, [15, 6, 7, 8, 9], first=d[15], last=d[9])
            _line_embedding(g, [9, 10, 19, 20, 21], first=d[9], last=d[21])
            _line_embedding(g, [45, 54, 53, 52, 51], first=d[45], last=d[51])
            _line_embedding(g, [51, 50, 49, 58, 57], first=d[51], last=d[57])
            _line_embedding(g, [51, 32, 33, 34, 61], first=d[51], last=d[61])
            _line_embedding(g, [61, 62, 41, 40, 27], first=d[61], last=d[27])
            _line_embedding(g, [9, 0, 1, 26, 27], first=d[9], last=d[27])
            _line_embedding(g, [27, 28, 47, 46, 45], first=d[27], last=d[45])
            _line_embedding(g, [15, 16, 17, 60, 61], first=d[15], last=d[61])

            # Top vertices
            _line_embedding(g, [2, 18, 42, 48], first=(-1, 7), last=(3, 7))

            return g

        elif embedding == 2:
            return g
        else:
            raise ValueError("The value of embedding must be 1 or 2.")

    def HouseGraph(self):
        """
        Returns a house graph with 5 nodes.

        A house graph is named for its shape. It is a triangle (roof) over a
        square (walls).

        This constructor depends on NetworkX numeric labeling.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the house
        graph is drawn with the first node in the lower-left corner of the
        house, the second in the lower-right corner of the house. The third
        node is in the upper-left corner connecting the roof to the wall,
        and the fourth is in the upper-right corner connecting the roof to
        the wall. The fifth node is the top of the roof, connected only to
        the third and fourth.

        EXAMPLES: Construct and show a house graph

        ::

            sage: g = graphs.HouseGraph()
            sage: g.show() # long time
        """
        pos_dict = {0:(-1,0),1:(1,0),2:(-1,1),3:(1,1),4:(0,2)}
        import networkx
        G = networkx.house_graph()
        return graph.Graph(G, pos=pos_dict, name="House Graph")

    def HouseXGraph(self):
        """
        Returns a house X graph with 5 nodes.

        A house X graph is a house graph with two additional edges. The
        upper-right corner is connected to the lower-left. And the
        upper-left corner is connected to the lower-right.

        This constructor depends on NetworkX numeric labeling.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the house X
        graph is drawn with the first node in the lower-left corner of the
        house, the second in the lower-right corner of the house. The third
        node is in the upper-left corner connecting the roof to the wall,
        and the fourth is in the upper-right corner connecting the roof to
        the wall. The fifth node is the top of the roof, connected only to
        the third and fourth.

        EXAMPLES: Construct and show a house X graph

        ::

            sage: g = graphs.HouseXGraph()
            sage: g.show() # long time
        """
        pos_dict = {0:(-1,0),1:(1,0),2:(-1,1),3:(1,1),4:(0,2)}
        import networkx
        G = networkx.house_x_graph()
        return graph.Graph(G, pos=pos_dict, name="House Graph")

    def KrackhardtKiteGraph(self):
        """
        Returns a Krackhardt kite graph with 10 nodes.

        The Krackhardt kite graph was originally developed by David
        Krackhardt for the purpose of studying social networks. It is used
        to show the distinction between: degree centrality, betweeness
        centrality, and closeness centrality. For more information read the
        plotting section below in conjunction with the example.

        REFERENCES:

        - [1] Kreps, V. (2002). "Social Network Analysis".  [Online] Available:
          http://www.fsu.edu/~spap/water/network/intro.htm [2007,
          January 17]

        This constructor depends on NetworkX numeric labeling.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the graph is
        drawn left to right, in top to bottom row sequence of [2, 3, 2, 1,
        1, 1] nodes on each row. This places the fourth node (3) in the
        center of the kite, with the highest degree. But the fourth node
        only connects nodes that are otherwise connected, or those in its
        clique (i.e.: Degree Centrality). The eighth (7) node is where the
        kite meets the tail. It has degree = 3, less than the average, but
        is the only connection between the kite and tail (i.e.: Betweenness
        Centrality). The sixth and seventh nodes (5 and 6) are drawn in the
        third row and have degree = 5. These nodes have the shortest path
        to all other nodes in the graph (i.e.: Closeness Centrality).
        Please execute the example for visualization.

        EXAMPLE: Construct and show a Krackhardt kite graph

        ::

            sage: g = graphs.KrackhardtKiteGraph()
            sage: g.show() # long time
        """
        pos_dict = {0:(-1,4),1:(1,4),2:(-2,3),3:(0,3),4:(2,3),5:(-1,2),6:(1,2),7:(0,1),8:(0,0),9:(0,-1)}
        import networkx
        G = networkx.krackhardt_kite_graph()
        return graph.Graph(G, pos=pos_dict, name="Krackhardt Kite Graph")

    def LadderGraph(self, n):
        """
        Returns a ladder graph with 2\*n nodes.

        A ladder graph is a basic structure that is typically displayed as
        a ladder, i.e.: two parallel path graphs connected at each
        corresponding node pair.

        This constructor depends on NetworkX numeric labels.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, each ladder
        graph will be displayed horizontally, with the first n nodes
        displayed left to right on the top horizontal line.

        EXAMPLES: Construct and show a ladder graph with 14 nodes

        ::

            sage: g = graphs.LadderGraph(7)
            sage: g.show() # long time

        Create several ladder graphs in a Sage graphics array

        ::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.LadderGraph(i+2)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time
        """
        pos_dict = {}
        for i in range(n):
            pos_dict[i] = (i,1)
        for i in range(n,2*n):
            x = i - n
            pos_dict[i] = (x,0)
        import networkx
        G = networkx.ladder_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Ladder graph")

    def LollipopGraph(self, n1, n2):
        """
        Returns a lollipop graph with n1+n2 nodes.

        A lollipop graph is a path graph (order n2) connected to a complete
        graph (order n1). (A barbell graph minus one of the bells).

        This constructor depends on NetworkX numeric labels.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the complete
        graph will be drawn in the lower-left corner with the (n1)th node
        at a 45 degree angle above the right horizontal center of the
        complete graph, leading directly into the path graph.

        EXAMPLES: Construct and show a lollipop graph Candy = 13, Stick =
        4

        ::

            sage: g = graphs.LollipopGraph(13,4)
            sage: g.show() # long time

        Create several lollipop graphs in a Sage graphics array

        ::

            sage: g = []
            sage: j = []
            sage: for i in range(6):
            ...    k = graphs.LollipopGraph(i+3,4)
            ...    g.append(k)
            ...
            sage: for i in range(2):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time
        """
        pos_dict = {}

        for i in range(n1):
            x = float(cos((pi/4) - ((2*pi)/n1)*i) - n2/2 - 1)
            y = float(sin((pi/4) - ((2*pi)/n1)*i) - n2/2 - 1)
            j = n1-1-i
            pos_dict[j] = (x,y)
        for i in range(n1, n1+n2):
            x = float(i - n1 - n2/2 + 1)
            y = float(i - n1 - n2/2 + 1)
            pos_dict[i] = (x,y)

        import networkx
        G = networkx.lollipop_graph(n1,n2)
        return graph.Graph(G, pos=pos_dict, name="Lollipop Graph")

    def PathGraph(self, n, pos=None):
        """
        Returns a path graph with n nodes. Pos argument takes a string
        which is either 'circle' or 'line', (otherwise the default is
        used). See the plotting section below for more detail.

        A path graph is a graph where all inner nodes are connected to
        their two neighbors and the two end-nodes are connected to their
        one inner neighbors. (i.e.: a cycle graph without the first and
        last node connected).

        This constructor depends on NetworkX numeric labels.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the graph may
        be drawn in one of two ways: The 'line' argument will draw the
        graph in a horizontal line (left to right) if there are less than
        11 nodes. Otherwise the 'line' argument will append horizontal
        lines of length 10 nodes below, alternating left to right and right
        to left. The 'circle' argument will cause the graph to be drawn in
        a cycle-shape, with the first node at the top and then about the
        circle in a clockwise manner. By default (without an appropriate
        string argument) the graph will be drawn as a 'circle' if 10 n 41
        and as a 'line' for all other n.

        EXAMPLES: Show default drawing by size: 'line': n 11

        ::

            sage: p = graphs.PathGraph(10)
            sage: p.show() # long time

        'circle': 10 n 41

        ::

            sage: q = graphs.PathGraph(25)
            sage: q.show() # long time

        'line': n 40

        ::

            sage: r = graphs.PathGraph(55)
            sage: r.show() # long time

        Override the default drawing::

            sage: s = graphs.PathGraph(5,'circle')
            sage: s.show() # long time
        """
        pos_dict = {}

        # Choose appropriate drawing pattern
        circle = False
        if pos == "circle": circle = True
        elif pos == "line": circle = False
        # Otherwise use default by size of n
        elif 10 < n < 41: circle = True

        # Draw 'circle'
        if circle:
            for i in range(n):
                x = float(cos((pi/2) + ((2*pi)/n)*i))
                y = float(sin((pi/2) + ((2*pi)/n)*i))
                pos_dict[i] = (x,y)
        # Draw 'line'
        else:
            counter = 0 # node index
            rem = n%10 # remainder to appear on last row
            rows = n//10 # number of rows (not counting last row)
            lr = True # left to right

            for i in range(rows): # note that rows doesn't include last row
                y = -i
                for j in range(10):
                    if lr:
                        x = j
                    else:
                        x = 9 - j
                    pos_dict[counter] = (x,y)
                    counter += 1
                if lr: lr = False
                else: lr = True
            y = -rows
            for j in range(rem): # last row
                if lr:
                    x = j
                else:
                    x = 9 - j
                pos_dict[counter] = (x,y)
                counter += 1

        import networkx
        G = networkx.path_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Path Graph")

    def StarGraph(self, n):
        """
        Returns a star graph with n+1 nodes.

        A Star graph is a basic structure where one node is connected to
        all other nodes.

        This constructor is dependent on NetworkX numeric labels.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, each star
        graph will be displayed with the first (0) node in the center, the
        second node (1) at the top, with the rest following in a
        counterclockwise manner. (0) is the node connected to all other
        nodes.

        The star graph is a good opportunity to compare efficiency of
        filling a position dictionary vs. using the spring-layout algorithm
        for plotting. As far as display, the spring-layout should push all
        other nodes away from the (0) node, and thus look very similar to
        this constructor's positioning.

        EXAMPLES::

            sage: import networkx

        Compare the plots::

            sage: n = networkx.star_graph(23)
            sage: spring23 = Graph(n)
            sage: posdict23 = graphs.StarGraph(23)
            sage: spring23.show() # long time
            sage: posdict23.show() # long time

        View many star graphs as a Sage Graphics Array

        With this constructor (i.e., the position dictionary filled)

        ::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.StarGraph(i+3)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time

        Compared to plotting with the spring-layout algorithm

        ::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    spr = networkx.star_graph(i+3)
            ...    k = Graph(spr)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time
        """
        pos_dict = {}
        pos_dict[0] = (0,0)
        for i in range(1,n+1):
            x = float(cos((pi/2) + ((2*pi)/n)*(i-1)))
            y = float(sin((pi/2) + ((2*pi)/n)*(i-1)))
            pos_dict[i] = (x,y)
        import networkx
        G = networkx.star_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Star graph")

    def WheelGraph(self, n):
        """
        Returns a Wheel graph with n nodes.

        A Wheel graph is a basic structure where one node is connected to
        all other nodes and those (outer) nodes are connected cyclically.

        This constructor depends on NetworkX numeric labels.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, each wheel
        graph will be displayed with the first (0) node in the center, the
        second node at the top, and the rest following in a
        counterclockwise manner.

        With the wheel graph, we see that it doesn't take a very large n at
        all for the spring-layout to give a counter-intuitive display. (See
        Graphics Array examples below).

        EXAMPLES: We view many wheel graphs with a Sage Graphics Array,
        first with this constructor (i.e., the position dictionary
        filled)::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.WheelGraph(i+3)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time

        Next, using the spring-layout algorithm::

            sage: import networkx
            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    spr = networkx.wheel_graph(i+3)
            ...    k = Graph(spr)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time

        Compare the plotting::

            sage: n = networkx.wheel_graph(23)
            sage: spring23 = Graph(n)
            sage: posdict23 = graphs.WheelGraph(23)
            sage: spring23.show() # long time
            sage: posdict23.show() # long time
        """
        pos_dict = {}
        pos_dict[0] = (0,0)
        for i in range(1,n):
            x = float(cos((pi/2) + ((2*pi)/(n-1))*(i-1)))
            y = float(sin((pi/2) + ((2*pi)/(n-1))*(i-1)))
            pos_dict[i] = (x,y)
        import networkx
        G = networkx.wheel_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Wheel graph")

################################################################################
#   Platonic Solids
################################################################################

    def TetrahedralGraph(self):
        """
        Returns a tetrahedral graph (with 4 nodes).

        A tetrahedron is a 4-sided triangular pyramid. The tetrahedral
        graph corresponds to the connectivity of the vertices of the
        tetrahedron. This graph is equivalent to a wheel graph with 4 nodes
        and also a complete graph on four nodes. (See examples below).

        PLOTTING: The tetrahedral graph should be viewed in 3 dimensions.
        We chose to use the default spring-layout algorithm here, so that
        multiple iterations might yield a different point of reference for
        the user. We hope to add rotatable, 3-dimensional viewing in the
        future. In such a case, a string argument will be added to select
        the flat spring-layout over a future implementation.

        EXAMPLES: Construct and show a Tetrahedral graph

        ::

            sage: g = graphs.TetrahedralGraph()
            sage: g.show() # long time

        The following example requires networkx::

            sage: import networkx as NX

        Compare this Tetrahedral, Wheel(4), Complete(4), and the
        Tetrahedral plotted with the spring-layout algorithm below in a
        Sage graphics array::

            sage: tetra_pos = graphs.TetrahedralGraph()
            sage: tetra_spring = Graph(NX.tetrahedral_graph())
            sage: wheel = graphs.WheelGraph(4)
            sage: complete = graphs.CompleteGraph(4)
            sage: g = [tetra_pos, tetra_spring, wheel, complete]
            sage: j = []
            sage: for i in range(2):
            ...    n = []
            ...    for m in range(2):
            ...        n.append(g[i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time
        """
        import networkx
        G = networkx.tetrahedral_graph()
        return graph.Graph(G, name="Tetrahedron", pos =
                           { 0 : (0, 0),
                             1 : (0, 1),
                             2 : (cos(3.5*pi/3), sin(3.5*pi/3)),
                             3 : (cos(5.5*pi/3), sin(5.5*pi/3))}
                           )

    def HexahedralGraph(self):
        """
        Returns a hexahedral graph (with 8 nodes).

        A regular hexahedron is a 6-sided cube. The hexahedral graph
        corresponds to the connectivity of the vertices of the hexahedron.
        This graph is equivalent to a 3-cube.

        PLOTTING: The hexahedral graph should be viewed in 3 dimensions. We
        chose to use the default spring-layout algorithm here, so that
        multiple iterations might yield a different point of reference for
        the user. We hope to add rotatable, 3-dimensional viewing in the
        future. In such a case, a string argument will be added to select
        the flat spring-layout over a future implementation.

        EXAMPLES: Construct and show a Hexahedral graph

        ::

            sage: g = graphs.HexahedralGraph()
            sage: g.show() # long time

        Create several hexahedral graphs in a Sage graphics array. They
        will be drawn differently due to the use of the spring-layout
        algorithm.

        ::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.HexahedralGraph()
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time
        """
        return graph.Graph({0:[1,3,4], 1:[2,5], 2:[3,6], 3:[7], 4:[5,7],\
                            5:[6], 6:[7]},
                           name="Hexahedron",
                           pos = {
                              0 : (0,0),
                              1 : (1,0),
                              3 : (0,1),
                              2 : (1,1),
                              4 : (.5,.5),
                              5 : (1.5,.5),
                              7 : (.5,1.5),
                              6 : (1.5,1.5)
                              })

    def OctahedralGraph(self):
        """
        Returns an Octahedral graph (with 6 nodes).

        The regular octahedron is an 8-sided polyhedron with triangular
        faces. The octahedral graph corresponds to the connectivity of the
        vertices of the octahedron. It is the line graph of the tetrahedral
        graph. The octahedral is symmetric, so the spring-layout algorithm
        will be very effective for display.

        PLOTTING: The Octahedral graph should be viewed in 3 dimensions. We
        chose to use the default spring-layout algorithm here, so that
        multiple iterations might yield a different point of reference for
        the user. We hope to add rotatable, 3-dimensional viewing in the
        future. In such a case, a string argument will be added to select
        the flat spring-layout over a future implementation.

        EXAMPLES: Construct and show an Octahedral graph

        ::

            sage: g = graphs.OctahedralGraph()
            sage: g.show() # long time

        Create several octahedral graphs in a Sage graphics array They will
        be drawn differently due to the use of the spring-layout algorithm

        ::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.OctahedralGraph()
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time
        """
        import networkx
        G = networkx.octahedral_graph()

        pos = {}
        r1 = 5
        r2 = 1
        for i,v in enumerate([0,1,2]):
            i = i + 0.75
            pos[v] = (r1*cos(i*2*pi/3),r1*sin(i*2*pi/3))

        for i,v in enumerate([4,3,5]):
            i = i + .25
            pos[v] = (r2*cos(i*2*pi/3),r2*sin(i*2*pi/3))


        return graph.Graph(G, name="Octahedron", pos=pos)

    def IcosahedralGraph(self):
        """
        Returns an Icosahedral graph (with 12 nodes).

        The regular icosahedron is a 20-sided triangular polyhedron. The
        icosahedral graph corresponds to the connectivity of the vertices
        of the icosahedron. It is dual to the dodecahedral graph. The
        icosahedron is symmetric, so the spring-layout algorithm will be
        very effective for display.

        PLOTTING: The Icosahedral graph should be viewed in 3 dimensions.
        We chose to use the default spring-layout algorithm here, so that
        multiple iterations might yield a different point of reference for
        the user. We hope to add rotatable, 3-dimensional viewing in the
        future. In such a case, a string argument will be added to select
        the flat spring-layout over a future implementation.

        EXAMPLES: Construct and show an Octahedral graph

        ::

            sage: g = graphs.IcosahedralGraph()
            sage: g.show() # long time

        Create several icosahedral graphs in a Sage graphics array. They
        will be drawn differently due to the use of the spring-layout
        algorithm.

        ::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.IcosahedralGraph()
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time
        """
        import networkx
        G = networkx.icosahedral_graph()

        pos = {}
        r1 = 5
        r2 = 2
        for i,v in enumerate([2,8,7,11,4,6]):
            i = i + .5
            pos[v] = (r1*cos(i*pi/3),r1*sin(i*pi/3))

        for i,v in enumerate([1,9,0,10,5,3]):
            i = i + .5
            pos[v] = (r2*cos(i*pi/3),r2*sin(i*pi/3))

        return graph.Graph(G, name="Icosahedron", pos = pos)

    def DodecahedralGraph(self):
        """
        Returns a Dodecahedral graph (with 20 nodes)

        The dodecahedral graph is cubic symmetric, so the spring-layout
        algorithm will be very effective for display. It is dual to the
        icosahedral graph.

        PLOTTING: The Dodecahedral graph should be viewed in 3 dimensions.
        We chose to use the default spring-layout algorithm here, so that
        multiple iterations might yield a different point of reference for
        the user. We hope to add rotatable, 3-dimensional viewing in the
        future. In such a case, a string argument will be added to select
        the flat spring-layout over a future implementation.

        EXAMPLES: Construct and show a Dodecahedral graph

        ::

            sage: g = graphs.DodecahedralGraph()
            sage: g.show() # long time

        Create several dodecahedral graphs in a Sage graphics array They
        will be drawn differently due to the use of the spring-layout
        algorithm

        ::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.DodecahedralGraph()
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time
        """
        import networkx
        G = networkx.dodecahedral_graph()

        pos = {}
        r1 = 7
        r2 = 4.7
        r3 = 3.8
        r4 = 1.5

        for i,v in enumerate([19,0,1,2,3]):
            i = i + .25
            pos[v] = (r1*cos(i*2*pi/5),r1*sin(i*2*pi/5))

        for i,v in enumerate([18,10,8,6,4]):
            i = i + .25
            pos[v] = (r2*cos(i*2*pi/5),r2*sin(i*2*pi/5))

        for i,v in enumerate([17,11,9,7,5]):
            i = i - .25
            pos[v] = (r3*cos(i*2*pi/5),r3*sin(i*2*pi/5))

        for i,v in enumerate([12,13,14,15,16]):
            i = i + .75
            pos[v] = (r4*cos(i*2*pi/5),r4*sin(i*2*pi/5))

        return graph.Graph(G, name="Dodecahedron", pos=pos)

    #######################################################################
    #   Named Graphs
    #######################################################################

    def Balaban10Cage(self, embedding=1):
        r"""
        Returns the Balaban 10-cage.

        The Balaban 10-cage is a 3-regular graph with 70 vertices and
        105 edges. See its :wikipedia:`Wikipedia page
        <Balaban_10-cage>`.

        The default embedding gives a deeper understanding of the
        graph's automorphism group. It is divided into 4 layers (each
        layer being a set of points at equal distance from the drawing's
        center). From outside to inside:

        - L1: The outer layer (vertices which are the furthest from the
          origin) is actually the disjoint union of two cycles of length
          10.

        - L2: The second layer is an independent set of 20 vertices.

        - L3: The third layer is a matching on 10 vertices.

        - L4: The inner layer (vertices which are the closest from the
          origin) is also the disjoint union of two cycles of length 10.

        This graph is not vertex-transitive, and its vertices are
        partitioned into 3 orbits: L2, L3, and the union of L1 of L4
        whose elements are equivalent.

        INPUT:

        - ``embedding`` -- two embeddings are available, and can be
          selected by setting ``embedding`` to be either 1 or 2.

        EXAMPLES::

            sage: g = graphs.Balaban10Cage()
            sage: g.girth()
            10
            sage: g.chromatic_number()
            2
            sage: g.diameter()
            6
            sage: g.is_hamiltonian()
            True
            sage: g.show(figsize=[10,10])   # long time

        TESTS::

            sage: graphs.Balaban10Cage(embedding='foo')
            Traceback (most recent call last):
            ...
            ValueError: The value of embedding must be 1 or 2.
        """

        L = [-9, -25, -19, 29, 13, 35, -13, -29, 19, 25, 9, -29, 29, 17, 33,
              21, 9,-13, -31, -9, 25, 17, 9, -31, 27, -9, 17, -19, -29, 27,
              -17, -9, -29, 33, -25,25, -21, 17, -17, 29, 35, -29, 17, -17,
              21, -25, 25, -33, 29, 9, 17, -27, 29, 19, -17, 9, -27, 31, -9,
              -17, -25, 9, 31, 13, -9, -21, -33, -17, -29, 29]

        g = graphs.LCFGraph(70, L, 1)
        g.name("Balaban 10-cage")

        if embedding == 2:
            return g
        elif embedding != 1:
            raise ValueError("The value of embedding must be 1 or 2.")

        L3 = [5, 24, 35, 46, 29, 40, 51, 34, 45, 56]
        _circle_embedding(g, L3, center=(0,0), radius = 4.3)

        L2  = [6, 4, 23, 25, 60, 36, 1, 47, 28, 30, 39, 41, 50, 52, 33, 9, 44,
                20, 55, 57]
        _circle_embedding(g, L2, center=(0,0), radius = 5, shift=-.5)


        L1a = [69, 68, 67, 66, 65, 64, 63, 62, 61, 0]
        L1b = [19, 18, 17, 16, 15, 14, 13, 12, 11, 10]
        _circle_embedding(g, L1a, center=(0,0), radius = 6, shift = 3.25)
        _circle_embedding(g, L1b, center=(0,0), radius = 6, shift = -1.25)

        L4a = [37, 2, 31, 38, 53, 32, 21, 54, 3, 22]
        _circle_embedding(g, L4a, center=(0,0), radius = 3, shift = 1.9)

        L4b = [26, 59, 48, 27, 42, 49, 8, 43, 58, 7]
        _circle_embedding(g, L4b, center=(0,0), radius = 3, shift = 1.1)

        return g

    def Balaban11Cage(self, embedding = 1):
        r"""
        Returns the Balaban 11-cage.

        For more information, see this :wikipedia:`Wikipedia article on
        the Balaban 11-cage <Balaban_11-cage>`.

        INPUT:

        - ``embedding`` -- three embeddings are available, and can be
          selected by setting ``embedding`` to be 1, 2, or 3.

          - The first embedding is the one appearing on page 9 of the
            Fifth Annual Graph Drawing Contest report [FAGDC]_. It
            separates vertices based on their eccentricity (see
            :meth:`eccentricity()
            <sage.graphs.generic_graph.GenericGraph.eccentricity>`).

          - The second embedding has been produced just for Sage and is
            meant to emphasize the automorphism group's 6 orbits.

          - The last embedding is the default one produced by the
            :meth:`LCFGraph` constructor.

        .. NOTE::

            The vertex labeling changes according to the value of
            ``embedding=1``.

        EXAMPLES:

        Basic properties::

            sage: g = graphs.Balaban11Cage()
            sage: g.order()
            112
            sage: g.size()
            168
            sage: g.girth()
            11
            sage: g.diameter()
            8
            sage: g.automorphism_group().cardinality()
            64

        Our many embeddings::

            sage: g1 = graphs.Balaban11Cage(embedding=1)
            sage: g2 = graphs.Balaban11Cage(embedding=2)
            sage: g3 = graphs.Balaban11Cage(embedding=3)
            sage: g1.show(figsize=[10,10])   # long time
            sage: g2.show(figsize=[10,10])   # long time
            sage: g3.show(figsize=[10,10])   # long time

        Proof that the embeddings are the same graph::

            sage: g1.is_isomorphic(g2) # g2 and g3 are obviously isomorphic
            True

        TESTS::

            sage: graphs.Balaban11Cage(embedding='xyzzy')
            Traceback (most recent call last):
            ...
            ValueError: The value of embedding must be 1, 2, or 3.

        REFERENCES:

        .. [FAGDC] Fifth Annual Graph Drawing Contest
           P. Eaded, J. Marks, P.Mutzel, S. North
           http://www.merl.com/papers/docs/TR98-16.pdf
        """
        if embedding == 1:
            pos_dict = {}
            for j in range(8):
                for i in range(8):
                    pos_dict[str(j) + str(i)]= [
                            0.8 * float(cos(2*((8*j + i)*pi/64 + pi/128))),
                            0.8 * float(sin(2*((8*j + i)*pi/64 + pi/128)))
                    ]
                for i in range(4):
                    pos_dict['1' + str(j) + str(i)] = [
                            1.1 * float(cos(2*((4*j + i)*pi/32 + pi/64))),
                            1.1 * float(sin(2*((4*j + i)*pi/32 + pi/64)))
                    ]
                for i in range(2):
                    pos_dict['1' + str(j) + str(i + 4)] = [
                            1.4 * float(cos(2*((2*j + i)*pi/16 + pi/32))),
                            1.4 * float(sin(2*((2*j + i)*pi/16 + pi/32)))
                    ]

            edge_dict = {
                "00": ["11"], "01": ["10"],   "02": ["53"], "03": ["52"],
                "11": ["20"], "10": ["21"],   "53": ["22"], "52": ["23"],
                "20": ["31"], "21": ["30"],   "22": ["33"], "23": ["32"],
                "31": ["40"], "30": ["41"],   "33": ["43"], "32": ["42"],
                "40": ["50"], "41": ["51"],   "43": ["12"], "42": ["13"],
                "50": ["61"], "51": ["60"],   "12": ["63"], "13": ["62"],
                "61": ["70"], "60": ["71"],   "63": ["72"], "62": ["73"],
                "70": ["01"], "71": ["00"],   "72": ["03"], "73": ["02"],

                "04": ["35"], "05": ["34"],   "06": ["37"], "07": ["36"],
                "35": ["64"], "34": ["65"],   "37": ["66"], "36": ["67"],
                "64": ["55"], "65": ["54"],   "66": ["17"], "67": ["16"],
                "55": ["45"], "54": ["44"],   "17": ["46"], "16": ["47"],
                "45": ["74"], "44": ["75"],   "46": ["76"], "47": ["77"],
                "74": ["25"], "75": ["24"],   "76": ["27"], "77": ["26"],
                "25": ["14"], "24": ["15"],   "27": ["56"], "26": ["57"],
                "14": ["05"], "15": ["04"],   "56": ["07"], "57": ["06"],

                "100": ["03", "04"],   "110": ["10", "12"],
                "101": ["01", "06"],   "111": ["11", "13"],
                "102": ["00", "07"],   "112": ["14", "16"],
                "103": ["02", "05"],   "113": ["15", "17"],

                "120": ["22", "24"],   "130": ["33", "36"],
                "121": ["20", "26"],   "131": ["32", "37"],
                "122": ["21", "27"],   "132": ["31", "34"],
                "123": ["23", "25"],   "133": ["30", "35"],

                "140": ["43", "45"],   "150": ["50", "52"],
                "141": ["40", "46"],   "151": ["51", "53"],
                "142": ["41", "47"],   "152": ["54", "56"],
                "143": ["42", "44"],   "153": ["55", "57"],

                "160": ["60", "66"],   "170": ["73", "76"],
                "161": ["63", "65"],   "171": ["72", "77"],
                "162": ["62", "64"],   "172": ["71", "74"],
                "163": ["61", "67"],   "173": ["70", "75"],

                "104": ["100", "102", "105"],   "114": ["110", "111", "115"],
                "105": ["101", "103", "104"],   "115": ["112", "113", "114"],

                "124": ["120", "121", "125"],   "134": ["130", "131", "135"],
                "125": ["122", "123", "124"],   "135": ["132", "133", "134"],

                "144": ["140", "141", "145"],   "154": ["150", "151", "155"],
                "145": ["142", "143", "144"],   "155": ["152", "153", "154"],

                "164": ["160", "161", "165"],   "174": ["170", "171", "175"],
                "165": ["162", "163", "164"],   "175": ["172", "173", "174"]
            }

            return graph.Graph(edge_dict, pos=pos_dict, name="Balaban 11-cage")

        elif embedding == 2 or embedding == 3:
            L = [44, 26, -47, -15, 35, -39, 11, -27, 38, -37, 43, 14, 28, 51,
                 -29, -16, 41, -11, -26, 15, 22, -51, -35, 36, 52, -14, -33,
                 -26, -46, 52, 26, 16, 43, 33, -15, 17, -53, 23, -42, -35, -28,
                 30, -22, 45, -44, 16, -38, -16, 50, -55, 20, 28, -17, -43,
                 47, 34, -26, -41, 11, -36, -23, -16, 41, 17, -51, 26, -33,
                 47, 17, -11, -20, -30, 21, 29, 36, -43, -52, 10, 39, -28, -17,
                 -52, 51, 26, 37, -17, 10, -10, -45, -34, 17, -26, 27, -21,
                 46, 53, -10, 29, -50, 35, 15, -47, -29, -41, 26, 33, 55, -17,
                 42, -26, -36, 16]

            g = graphs.LCFGraph(112, L, 1)
            g.name("Balaban 11-cage")

            if embedding == 3:
                return g

            v1 = [34, 2, 54, 43, 66, 20, 89, 100, 72, 76, 6, 58, 16, 78, 74,
                  70, 36, 94, 27, 25, 10, 8, 45, 60, 14, 64, 80, 82, 109, 107,
                  49, 98]
            v2 = [88, 3, 19, 55, 67, 42, 101, 33, 77, 5, 17, 57, 69, 71, 73,
                  75, 11, 61, 28, 9, 37, 26, 46, 95, 13, 63, 81, 83, 108, 106,
                  48, 97]
            l1 = [35, 93, 1, 24, 53, 7, 44, 59, 15, 65, 79, 21, 110, 90, 50,
                  99]
            l2 = [87, 4, 18, 56, 68, 41, 102, 32, 12, 62, 29, 84, 38, 105, 47,
                  96]

            d = g.get_pos()
            for i,v in enumerate(v1):
                d[v] = (-2, 16.5-i)

            for i,v in enumerate(l1):
                d[v] = (-10, 8-i)

            for i,v in enumerate(l2):
                d[v] = (10, 8.5-i)

            for i,v in enumerate(v2):
                d[v] = (2, 16.5-i)

            for i,v in enumerate([0, 111, 92, 91, 52, 51, 23, 22]):
                d[v] = (-20, 14.5-4*i)

            for i,v in enumerate([104, 103, 86, 85, 40, 39, 31, 30]):
                d[v] = (20, 14.5-4*i)

            return g

        else:
            raise ValueError("The value of embedding must be 1, 2, or 3.")

    def BidiakisCube(self):
        r"""
        Returns the Bidiakis cube.

        For more information, see this
        `Wikipedia article on the Bidiakis cube <http://en.wikipedia.org/wiki/Bidiakis_cube>`_.

        EXAMPLES:

        The Bidiakis cube is a 3-regular graph having 12 vertices and 18
        edges. This means that each vertex has a degree of 3. ::

            sage: g = graphs.BidiakisCube(); g
            Bidiakis cube: Graph on 12 vertices
            sage: g.show()  # long time
            sage: g.order()
            12
            sage: g.size()
            18
            sage: g.is_regular(3)
            True

        It is a Hamiltonian graph with diameter 3 and girth 4::

            sage: g.is_hamiltonian()
            True
            sage: g.diameter()
            3
            sage: g.girth()
            4

        It is a planar graph with characteristic polynomial
        `(x - 3) (x - 2) (x^4) (x + 1) (x + 2) (x^2 + x - 4)^2` and
        chromatic number 3::

            sage: g.is_planar()
            True
            sage: bool(g.characteristic_polynomial() == expand((x - 3) * (x - 2) * (x^4) * (x + 1) * (x + 2) * (x^2 + x - 4)^2))
            True
            sage: g.chromatic_number()
            3
        """
        edge_dict = {
            0:[1,6,11], 1:[2,5], 2:[3,10], 3:[4,9], 4:[5,8],
            5:[6], 6:[7], 7:[8,11], 8:[9], 9:[10], 10:[11]}
        pos_dict = {
            0: [0, 1],
            1: [0.5, 0.866025403784439],
            2: [0.866025403784439, 0.500000000000000],
            3: [1, 0],
            4: [0.866025403784439, -0.5],
            5: [0.5, -0.866025403784439],
            6: [0, -1],
            7: [-0.5, -0.866025403784439],
            8: [-0.866025403784439, -0.5],
            9: [-1, 0],
            10: [-0.866025403784439, 0.5],
            11: [-0.5, 0.866025403784439]}
        return graph.Graph(edge_dict, pos=pos_dict, name="Bidiakis cube")

    def BiggsSmithGraph(self, embedding=1):
        r"""
        Returns the Biggs-Smith graph.

        For more information, see this :wikipedia:`Wikipedia article on
        the Biggs-Smith graph <Biggs-Smith_graph>`.

        INPUT:

        - ``embedding`` -- two embeddings are available, and can be
          selected by setting ``embedding`` to be 1 or 2.

        EXAMPLES:

        Basic properties::

            sage: g = graphs.BiggsSmithGraph()
            sage: g.order()
            102
            sage: g.size()
            153
            sage: g.girth()
            9
            sage: g.diameter()
            7
            sage: g.automorphism_group().cardinality()
            2448
            sage: g.show(figsize=[10, 10])   # long time

        The other embedding::

            sage: graphs.BiggsSmithGraph(embedding=2).show()

        TESTS::

            sage: graphs.BiggsSmithGraph(embedding='xyzzy')
            Traceback (most recent call last):
            ...
            ValueError: The value of embedding must be 1 or 2.

        """
        L = [16, 24, -38, 17, 34, 48, -19, 41, -35, 47, -20, 34, -36,
             21, 14, 48, -16, -36, -43, 28, -17, 21, 29, -43, 46, -24,
             28, -38, -14, -50, -45, 21, 8, 27, -21, 20, -37, 39, -34,
             -44, -8, 38, -21, 25, 15, -34, 18, -28, -41, 36, 8, -29,
             -21, -48, -28, -20, -47, 14, -8, -15, -27, 38, 24, -48, -18,
             25, 38, 31, -25, 24, -46, -14, 28, 11, 21, 35, -39, 43, 36,
             -38, 14, 50, 43, 36, -11, -36, -24, 45, 8, 19, -25, 38, 20,
             -24, -14, -21, -8, 44, -31, -38, -28, 37]

        g = graphs.LCFGraph(102, L, 1)
        g.name("Biggs-Smith graph")

        if embedding == 1:

            orbs = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 0],
                    [17, 101, 25, 66, 20, 38, 53, 89, 48, 75, 56, 92, 45, 78,
                     34, 28, 63],
                    [18, 36, 26, 65, 19, 37, 54, 90, 47, 76, 55, 91, 46, 77,
                     35, 27, 64],
                    [21, 39, 52, 88, 49, 74, 57, 93, 44, 79, 33, 29, 62, 83,
                     100, 24, 67],
                    [22, 97, 51, 96, 50, 95, 58, 94, 59, 80, 60, 81, 61, 82,
                     99, 23, 98],
                    [30, 86, 84, 72, 70, 68, 42, 40, 31, 87, 85, 73, 71, 69,
                     43, 41, 32]]

            # central orbits
            _circle_embedding(g, orbs[1], center=(-.4, 0), radius=.2)
            _circle_embedding(g, orbs[3], center=(.4, 0), radius=.2, shift=4)

            # lower orbits
            _circle_embedding(g, orbs[0], center=(-.9, -.5), radius=.3,
                    shift=2)
            _circle_embedding(g, orbs[2], center=(-.9, .5), radius=.3)

            # upper orbits
            _circle_embedding(g, orbs[4], center=(.9, -.5), radius=.3, shift=4)
            _circle_embedding(g, orbs[5], center=(.9, .5), radius=.3, shift=-2)

        elif embedding == 2:
            pass
        else:
            raise ValueError("The value of embedding must be 1 or 2.")

        return g

    def BrinkmannGraph(self):
        r"""
        Returns the Brinkmann graph.

        For more information, see the
        `Wikipedia article on the Brinkmann graph <http://en.wikipedia.org/wiki/Brinkmann_graph>`_.

        EXAMPLES:

        The Brinkmann graph is a 4-regular graph having 21 vertices and 42
        edges. This means that each vertex has degree 4. ::

            sage: G = graphs.BrinkmannGraph(); G
            Brinkmann graph: Graph on 21 vertices
            sage: G.show()  # long time
            sage: G.order()
            21
            sage: G.size()
            42
            sage: G.is_regular(4)
            True

        It is an Eulerian graph with radius 3, diameter 3, and girth 5. ::

            sage: G.is_eulerian()
            True
            sage: G.radius()
            3
            sage: G.diameter()
            3
            sage: G.girth()
            5

        The Brinkmann graph is also Hamiltonian with chromatic number 4::

            sage: G.is_hamiltonian()
            True
            sage: G.chromatic_number()
            4

        Its automorphism group is isomorphic to `D_7`::

            sage: ag = G.automorphism_group()
            sage: ag.is_isomorphic(DihedralGroup(7))
            True
        """
        edge_dict = {
            0: [2,5,7,13],
            1: [3,6,7,8],
            2: [4,8,9],
            3: [5,9,10],
            4: [6,10,11],
            5: [11,12],
            6: [12,13],
            7: [15,20],
            8: [14,16],
            9: [15,17],
            10: [16,18],
            11: [17,19],
            12: [18,20],
            13: [14,19],
            14: [17,18],
            15: [18,19],
            16: [19,20],
            17: [20]}
        pos_dict = {
            0: [0, 4],
            1: [3.12732592987212, 2.49395920743493],
            2: [3.89971164872729, -0.890083735825258],
            3: [1.73553495647023, -3.60387547160968],
            4: [-1.73553495647023, -3.60387547160968],
            5: [-3.89971164872729, -0.890083735825258],
            6: [-3.12732592987212, 2.49395920743493],
            7: [0.867767478235116, 1.80193773580484],
            8: [1.94985582436365, 0.445041867912629],
            9: [1.56366296493606, -1.24697960371747],
            10: [0, -2],
            11: [-1.56366296493606, -1.24697960371747],
            12: [-1.94985582436365, 0.445041867912629],
            13: [-0.867767478235116, 1.80193773580484],
            14: [0.433883739117558, 0.900968867902419],
            15: [0.974927912181824, 0.222520933956314],
            16: [0.781831482468030, -0.623489801858733],
            17: [0, -1],
            18: [-0.781831482468030, -0.623489801858733],
            19: [-0.974927912181824, 0.222520933956315],
            20: [-0.433883739117558, 0.900968867902419]}
        return graph.Graph(edge_dict, pos=pos_dict, name="Brinkmann graph")

    def DoubleStarSnark(self):
        r"""
        Returns the double star snark.

        The double star snark is a 3-regular graph on 30 vertices. See
        the :wikipedia:`Wikipedia page on the double star snark
        <Double-star_snark>`.

        EXAMPLES::

            sage: g = graphs.DoubleStarSnark()
            sage: g.order()
            30
            sage: g.size()
            45
            sage: g.chromatic_number()
            3
            sage: g.is_hamiltonian()
            False
            sage: g.automorphism_group().cardinality()
            80
            sage: g.show()
        """

        d = { 0: [1, 14, 15]
            , 1: [0, 2, 11]
            , 2: [1, 3, 7]
            , 3: [2, 4, 18]
            , 4: [3, 5, 14]
            , 5: [10, 4, 6]
            , 6: [5, 21, 7]
            , 7: [8, 2, 6]
            , 8: [9, 13, 7]
            , 9: [24, 8, 10]
            , 10: [9, 11, 5]
            , 11: [1, 10, 12]
            , 12: [11, 27, 13]
            , 13: [8, 12, 14]
            , 14: [0, 4, 13]
            , 15: [0, 16, 29]
            , 16: [15, 20, 23]
            , 17: [25, 18, 28]
            , 18: [3, 17, 19]
            , 19: [18, 26, 23]
            , 20: [16, 28, 21]
            , 21: [20, 6, 22]
            , 22: [26, 21, 29]
            , 23: [16, 24, 19]
            , 24: [25, 9, 23]
            , 25: [24, 17, 29]
            , 26: [27, 19, 22]
            , 27: [12, 26, 28]
            , 28: [17, 27, 20]
            , 29: [25, 22, 15]
            }

        g = graph.Graph(d, pos={}, name="Double star snark")
        _circle_embedding(g, range(15), radius=2)
        _circle_embedding(g, range(15, 30), radius=1.4)

        return g

    def ChvatalGraph(self):
        r"""
        Returns the Chvatal graph.

        Chvatal graph is one of the few known graphs to satisfy Grunbaum's
        conjecture that for every m, n, there is an m-regular,
        m-chromatic graph of girth at least n. For more information, see this
        `Wikipedia article on the Chvatal graph <http://en.wikipedia.org/wiki/Chv%C3%A1tal_graph>`_.

        EXAMPLES:

        The Chvatal graph has 12 vertices and 24 edges. It is a 4-regular,
        4-chromatic graph with radius 2, diameter 2, and girth 4. ::

            sage: G = graphs.ChvatalGraph(); G
            Chvatal graph: Graph on 12 vertices
            sage: G.order(); G.size()
            12
            24
            sage: G.degree()
            [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
            sage: G.chromatic_number()
            4
            sage: G.radius(); G.diameter(); G.girth()
            2
            2
            4
        """
        import networkx
        pos_dict = {}
        for i in range(5, 10):
            x = float(cos((pi / 2) + ((2 * pi) / 5) * i))
            y = float(sin((pi / 2) + ((2 * pi) / 5) * i))
            pos_dict[i] = (x, y)
        for i in range(5):
            x = float(2 * (cos((pi / 2) + ((2 * pi) / 5) * (i - 5))))
            y = float(2 * (sin((pi / 2) + ((2 * pi) / 5) * (i - 5))))
            pos_dict[i] = (x, y)
        pos_dict[10] = (0.5, 0)
        pos_dict[11] = (-0.5, 0)

        return graph.Graph(networkx.chvatal_graph(), pos=pos_dict, name="Chvatal graph")

    def ClebschGraph(self):
        r"""
        Return the Clebsch graph.

        EXAMPLES::

            sage: g = graphs.ClebschGraph()
            sage: g.automorphism_group().cardinality()
            1920
            sage: g.girth()
            4
            sage: g.chromatic_number()
            4
            sage: g.diameter()
            2
            sage: g.show(figsize=[10, 10]) # long time
        """
        g = graph.Graph(pos={})
        x = 0
        for i in range(8):
            g.add_edge(x % 16, (x + 1) % 16)
            g.add_edge(x % 16, (x + 6) % 16)
            g.add_edge(x % 16, (x + 8) % 16)
            x += 1
            g.add_edge(x % 16, (x + 3) % 16)
            g.add_edge(x % 16, (x + 2) % 16)
            g.add_edge(x % 16, (x + 8) % 16)
            x += 1

        _circle_embedding(g, range(16), shift=.5)
        g.name("Clebsch graph")

        return g

    def CoxeterGraph(self):
        r"""
        Return the Coxeter graph.

        See the :wikipedia:`Wikipedia page on the Coxeter graph
        <Coxeter_graph>`.

        EXAMPLES::

            sage: g = graphs.CoxeterGraph()
            sage: g.automorphism_group().cardinality()
            336
            sage: g.girth()
            7
            sage: g.chromatic_number()
            3
            sage: g.diameter()
            4
            sage: g.show(figsize=[10, 10]) # long time
        """
        g = graph.Graph({
                27: [6, 22, 14],
                24: [0, 7, 18],
                25: [8, 15, 2],
                26: [10, 16, 23],
                }, pos={})

        g.add_cycle(range(24))
        g.add_edges([(5, 11), (9, 20), (12, 1), (13, 19), (17, 4), (3, 21)])

        _circle_embedding(g, range(24))
        _circle_embedding(g, [24, 25, 26], radius=.5)
        g.get_pos()[27] = (0, 0)

        g.name("Coxeter Graph")

        return g

    def DesarguesGraph(self):
        """
        Returns the Desargues graph.

        PLOTTING: The layout chosen is the same as on the cover of [1].

        EXAMPLE::

            sage: D = graphs.DesarguesGraph()
            sage: L = graphs.LCFGraph(20,[5,-5,9,-9],5)
            sage: D.is_isomorphic(L)
            True
            sage: D.show()  # long time

        REFERENCE:

        - [1] Harary, F. Graph Theory. Reading, MA: Addison-Wesley,
          1994.
        """
        G=graphs.GeneralizedPetersenGraph(10,3)
        G.name("Desargues Graph")
        return G

    def DurerGraph(self):
        r"""
        Returns the Dürer graph.

        For more information, see this
        `Wikipedia article on the Dürer graph <http://en.wikipedia.org/wiki/D%C3%BCrer_graph>`_.

        EXAMPLES:

        The Dürer graph is named after Albrecht Dürer. It is a planar graph
        with 12 vertices and 18 edges. ::

            sage: G = graphs.DurerGraph(); G
            Durer graph: Graph on 12 vertices
            sage: G.is_planar()
            True
            sage: G.order()
            12
            sage: G.size()
            18

        The Dürer graph has chromatic number 3, diameter 4, and girth 3. ::

            sage: G.chromatic_number()
            3
            sage: G.diameter()
            4
            sage: G.girth()
            3

        Its automorphism group is isomorphic to `D_6`. ::

            sage: ag = G.automorphism_group()
            sage: ag.is_isomorphic(DihedralGroup(6))
            True
        """
        edge_dict = {
            0: [1,5,6],
            1: [2,7],
            2: [3,8],
            3: [4,9],
            4: [5,10],
            5: [11],
            6: [8,10],
            7: [9,11],
            8: [10],
            9: [11]}
        pos_dict = {
            0: [2, 0],
            1: [1, 1.73205080756888],
            2: [-1, 1.73205080756888],
            3: [-2, 0],
            4: [-1, -1.73205080756888],
            5: [1, -1.73205080756888],
            6: [1, 0],
            7: [0.5, 0.866025403784439],
            8: [-0.5, 0.866025403784439],
            9: [-1, 0],
            10: [-0.5, -0.866025403784439],
            11: [0.5, -0.866025403784439]}
        return graph.Graph(edge_dict, pos=pos_dict, name="Durer graph")

    def DyckGraph(self):
        """
        Returns the Dyck graph.

        For more information, see the `MathWorld article on the Dyck graph
        <http://mathworld.wolfram.com/DyckGraph.html>`_ or the `Wikipedia
        article on the Dyck graph <http://en.wikipedia.org/wiki/Dyck_graph>`_.

        EXAMPLES:

        The Dyck graph was defined by Walther von Dyck in 1881. It has `32`
        vertices and `48` edges, and is a cubic graph (regular of degree `3`)::

            sage: G = graphs.DyckGraph(); G
            Dyck graph: Graph on 32 vertices
            sage: G.order()
            32
            sage: G.size()
            48
            sage: G.is_regular()
            True
            sage: G.is_regular(3)
            True

        It is non-planar and Hamiltonian, as well as bipartite (making it a
        bicubic graph)::

            sage: G.is_planar()
            False
            sage: G.is_hamiltonian()
            True
            sage: G.is_bipartite()
            True

        It has radius `5`, diameter `5`, and girth `6`::

            sage: G.radius()
            5
            sage: G.diameter()
            5
            sage: G.girth()
            6

        Its chromatic number is `2` and its automorphism group is of order
        `192`::

            sage: G.chromatic_number()
            2
            sage: G.automorphism_group().cardinality()
            192

        It is a non-integral graph as it has irrational eigenvalues::

            sage: G.characteristic_polynomial().factor()
            (x - 3) * (x + 3) * (x - 1)^9 * (x + 1)^9 * (x^2 - 5)^6

        It is a toroidal graph, and its embedding on a torus is dual to an
        embedding of the Shrikhande graph (:meth:`ShrikhandeGraph
        <GraphGenerators.ShrikhandeGraph>`).
        """
        pos_dict = {}
        for i in range(8):
            pos_dict[i] = [float(cos((2*i) * pi/8)),
                           float(sin((2*i) * pi/8))]
            pos_dict[8 + i]  = [0.75 * pos_dict[i][0],
                                0.75 * pos_dict[i][1]]
            pos_dict[16 + i] = [0.50 * pos_dict[i][0],
                                0.50 * pos_dict[i][1]]
            pos_dict[24 + i] = [0.25 * pos_dict[i][0],
                                0.25 * pos_dict[i][1]]

        edge_dict = {
            0O00: [0O07, 0O01,   0O10], 0O10: [0O00,   0O27, 0O21],
            0O01: [0O00, 0O02,   0O11], 0O11: [0O01,   0O20, 0O22],
            0O02: [0O01, 0O03,   0O12], 0O12: [0O02,   0O21, 0O23],
            0O03: [0O02, 0O04,   0O13], 0O13: [0O03,   0O22, 0O24],
            0O04: [0O03, 0O05,   0O14], 0O14: [0O04,   0O23, 0O25],
            0O05: [0O04, 0O06,   0O15], 0O15: [0O05,   0O24, 0O26],
            0O06: [0O05, 0O07,   0O16], 0O16: [0O06,   0O25, 0O27],
            0O07: [0O06, 0O00,   0O17], 0O17: [0O07,   0O26, 0O20],

            0O20: [0O17, 0O11,   0O30], 0O30: [0O20,   0O35, 0O33],
            0O21: [0O10, 0O12,   0O31], 0O31: [0O21,   0O36, 0O34],
            0O22: [0O11, 0O13,   0O32], 0O32: [0O22,   0O37, 0O35],
            0O23: [0O12, 0O14,   0O33], 0O33: [0O23,   0O30, 0O36],
            0O24: [0O13, 0O15,   0O34], 0O34: [0O24,   0O31, 0O37],
            0O25: [0O14, 0O16,   0O35], 0O35: [0O25,   0O32, 0O30],
            0O26: [0O15, 0O17,   0O36], 0O36: [0O26,   0O33, 0O31],
            0O27: [0O16, 0O10,   0O37], 0O37: [0O27,   0O34, 0O32],
        }

        return graph.Graph(edge_dict, pos=pos_dict, name="Dyck graph")

    def EllinghamHorton54Graph(self):
        r"""
        Returns the Ellingham-Horton 54-graph.

        For more information, see the :wikipedia:`Wikipedia page on the
        Ellingham-Horton graphs <Ellingham-Horton_graph>`

        EXAMPLE:

        This graph is 3-regular::

            sage: g = graphs.EllinghamHorton54Graph()
            sage: g.is_regular(k=3)
            True

        It is 3-connected and bipartite::

            sage: g.vertex_connectivity() # not tested - too long
            3
            sage: g.is_bipartite()
            True

        It is not Hamiltonian::

            sage: g.is_hamiltonian() # not tested - too long
            False

        ... and it has a nice drawing ::

            sage: g.show(figsize=[10, 10]) # not tested - too long

        TESTS::

            sage: g.show() # long time
        """
        up = graphs.CycleGraph(16)
        low = 2*graphs.CycleGraph(6)

        for v in range(6):
            low.add_edge(v, v + 12)
            low.add_edge(v + 6, v + 12)
        low.add_edge(12, 15)
        low.delete_edge(1, 2)
        low.delete_edge(8, 7)
        low.add_edge(1, 8)
        low.add_edge(7, 2)


        # The set of vertices on top is 0..15
        # Bottom left is 16..33
        # Bottom right is 34..52
        # The two other vertices are 53, 54
        g = up + 2*low
        g.name("Ellingham-Horton 54-graph")
        g.set_pos({})

        g.add_edges([(15, 4), (3, 8), (7, 12), (11, 0), (2, 13), (5, 10)])
        g.add_edges([(30, 6), (29, 9), (48, 14), (47, 1)])
        g.add_edge(32, 52)
        g.add_edge(50, 52)
        g.add_edge(33, 53)
        g.add_edge(51, 53)
        g.add_edge(52, 53)

        # Top
        _circle_embedding(g, range(16), center=(0, .5), shift=.5, radius=.5)

        # Bottom-left
        _circle_embedding(g, range(16, 22), center=(-1.5, -1))
        _circle_embedding(g, range(22, 28), center=(-1.5, -1), radius=.5)
        _circle_embedding(g, range(28, 34), center=(-1.5, -1), radius=.7)

        # Bottom right
        _circle_embedding(g, range(34, 40), center=(1.5, -1))
        _circle_embedding(g, range(40, 46), center=(1.5, -1), radius=.5)
        _circle_embedding(g, range(46, 52), center=(1.5, -1), radius=.7)

        d = g.get_pos()
        d[52] = (-.3, -2.5)
        d[53] = (.3, -2.5)
        d[31] = (-2.2, -.9)
        d[28] = (-.8, -.9)
        d[46] = (2.2, -.9)
        d[49] = (.8, -.9)


        return g

    def EllinghamHorton78Graph(self):
        r"""
        Returns the Ellingham-Horton 78-graph.

        For more information, see the :wikipedia:`Wikipedia page on the
        Ellingham-Horton graphs
        <http://en.wikipedia.org/wiki/Ellingham%E2%80%93Horton_graph>`

        EXAMPLE:

        This graph is 3-regular::

            sage: g = graphs.EllinghamHorton78Graph()
            sage: g.is_regular(k=3)
            True

        It is 3-connected and bipartite::

            sage: g.vertex_connectivity() # not tested - too long
            3
            sage: g.is_bipartite()
            True

        It is not Hamiltonian::

            sage: g.is_hamiltonian() # not tested - too long
            False

        ... and it has a nice drawing ::

            sage: g.show(figsize=[10,10]) # not tested - too long

        TESTS::

            sage: g.show(figsize=[10, 10]) # not tested - too long
        """
        g = graph.Graph({
                0: [1, 5, 60], 1: [2, 12], 2: [3, 7], 3: [4, 14], 4: [5, 9],
                5: [6], 6: [7, 11], 7: [15], 8: [9, 13, 22], 9: [10],
                10: [11, 72], 11: [12], 12: [13], 13: [14], 14: [72],
                15: [16, 20], 16: [17, 27], 17: [18, 22], 18: [19, 29],
                19: [20, 24], 20: [21], 21: [22, 26], 23: [24, 28, 72],
                24: [25], 25: [26, 71], 26: [27], 27: [28], 28: [29],
                29: [69], 30: [31, 35, 52], 31: [32, 42], 32: [33, 37],
                33: [34, 43], 34: [35, 39], 35: [36], 36: [41, 63],
                37: [65, 66], 38: [39, 59, 74], 39: [40], 40: [41, 44],
                41: [42], 42: [74], 43: [44, 74], 44: [45], 45: [46, 50],
                46: [47, 57], 47: [48, 52], 48: [49, 75], 49: [50, 54],
                50: [51], 51: [52, 56], 53: [54, 58, 73], 54: [55],
                55: [56, 59], 56: [57], 57: [58], 58: [75], 59: [75],
                60: [61, 64], 61: [62, 71], 62: [63, 77], 63: [67],
                64: [65, 69], 65: [77], 66: [70, 73], 67: [68, 73],
                68: [69, 76], 70: [71, 76], 76: [77]}, pos={})

        _circle_embedding(g, range(15), center=(-2.5, 1.5))
        _circle_embedding(g, range(15, 30), center=(-2.5, -1.5))
        _circle_embedding(g, [30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
            42, 74, 43, 44], center=(2.5, 1.5))
        _circle_embedding(g, [45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56,
            57, 58, 75, 59], center=(2.5, -1.5))

        d = g.get_pos()

        d[76] = (-.2, -.1)
        d[77] = (.2, .1)
        d[38] = (2.2, .1)
        d[52] = (2.3, -.1)
        d[15] = (-2.1, -.1)
        d[72] = (-2.1, .1)

        _line_embedding(g, [60, 61, 62, 63], first=(-1, 2), last=(1, 2))
        _line_embedding(g, [64, 65, 37], first=(-.5, 1.5), last=(1.2, 1.5))
        _line_embedding(g, [66, 73, 67, 68, 69], first=(1.2, -2),
                last=(-.8, -2))
        _line_embedding(g, [66, 70, 71], first=(.7, -1.5), last=(-1, -1.5))

        g.name("Ellingham-Horton 78-graph")

        return g

    def ErreraGraph(self):
        r"""
        Returns the Errera graph.

        For more information, see this
        `Wikipedia article on the Errera graph <http://en.wikipedia.org/wiki/Errera_graph>`_.

        EXAMPLES:

        The Errera graph is named after Alfred Errera. It is a planar graph
        on 17 vertices and having 45 edges. ::

            sage: G = graphs.ErreraGraph(); G
            Errera graph: Graph on 17 vertices
            sage: G.is_planar()
            True
            sage: G.order()
            17
            sage: G.size()
            45

        The Errera graph is Hamiltonian with radius 3, diameter 4, girth 3,
        and chromatic number 4. ::

            sage: G.is_hamiltonian()
            True
            sage: G.radius()
            3
            sage: G.diameter()
            4
            sage: G.girth()
            3
            sage: G.chromatic_number()
            4

        Each vertex degree is either 5 or 6. That is, if `f` counts the
        number of vertices of degree 5 and `s` counts the number of vertices
        of degree 6, then `f + s` is equal to the order of the Errera
        graph. ::

            sage: D = G.degree_sequence()
            sage: D.count(5) + D.count(6) == G.order()
            True

        The automorphism group of the Errera graph is isomorphic to the
        dihedral group of order 20. ::

            sage: ag = G.automorphism_group()
            sage: ag.is_isomorphic(DihedralGroup(10))
            True
        """
        edge_dict = {
            0: [1,7,14,15,16],
            1: [2,9,14,15],
            2: [3,8,9,10,14],
            3: [4,9,10,11],
            4: [5,10,11,12],
            5: [6,11,12,13],
            6: [7,8,12,13,16],
            7: [13,15,16],
            8: [10,12,14,16],
            9: [11,13,15],
            10: [12],
            11: [13],
            13: [15],
            14: [16]}
        return graph.Graph(edge_dict, name="Errera graph")

    def FlowerSnark(self):
        """
        Returns a Flower Snark.

        A flower snark has 20 vertices. It is part of the class of
        biconnected cubic graphs with edge chromatic number = 4, known as
        snarks. (i.e.: the Petersen graph). All snarks are not Hamiltonian,
        non-planar and have Petersen graph graph minors.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the nodes are
        drawn 0-14 on the outer circle, and 15-19 in an inner pentagon.

        REFERENCES:

        - [1] Weisstein, E. (1999). "Flower Snark - from Wolfram
          MathWorld". [Online] Available:
          http://mathworld.wolfram.com/FlowerSnark.html [2007, February 17]

        EXAMPLES: Inspect a flower snark::

            sage: F = graphs.FlowerSnark()
            sage: F
            Flower Snark: Graph on 20 vertices
            sage: F.graph6_string()
            'ShCGHC@?GGg@?@?Gp?K??C?CA?G?_G?Cc'

        Now show it::

            sage: F.show() # long time
        """
        pos_dict = {}
        for i in range(15):
            x = float(2.5*(cos((pi/2) + ((2*pi)/15)*i)))
            y = float(2.5*(sin((pi/2) + ((2*pi)/15)*i)))
            pos_dict[i] = (x,y)
        for i in range(15,20):
            x = float(cos((pi/2) + ((2*pi)/5)*i))
            y = float(sin((pi/2) + ((2*pi)/5)*i))
            pos_dict[i] = (x,y)
        return graph.Graph({0:[1,14,15],1:[2,11],2:[3,7],3:[2,4,16],4:[5,14], \
                            5:[6,10],6:[5,7,17],8:[7,9,13],9:[10,18],11:[10,12], \
                            12:[13,19],13:[14],15:[19],16:[15,17],18:[17,19]}, \
                            pos=pos_dict, name="Flower Snark")

    def FosterGraph(self):
        """
        Returns the Foster graph.

        See the :wikipedia:`Wikipedia page on the Foster Graph
        <Foster_graph>`.

        EXAMPLE::

            sage: g = graphs.FosterGraph()
            sage: g.order()
            90
            sage: g.size()
            135
            sage: g.diameter()
            8
            sage: g.girth()
            10
            sage: g.automorphism_group().cardinality()
            4320
            sage: g.is_hamiltonian()
            True
        """
        g= graphs.LCFGraph(90, [17, -9, 37, -37, 9, -17], 15)
        g.name("Foster Graph")
        return g


    def FranklinGraph(self):
        r"""
        Returns the Franklin graph.

        For more information, see this
        `Wikipedia article on the Franklin graph <http://en.wikipedia.org/wiki/Franklin_graph>`_.

        EXAMPLES:

        The Franklin graph is named after Philip Franklin. It is a
        3-regular graph on 12 vertices and having 18 edges. ::

            sage: G = graphs.FranklinGraph(); G
            Franklin graph: Graph on 12 vertices
            sage: G.is_regular(3)
            True
            sage: G.order()
            12
            sage: G.size()
            18

        The Franklin graph is a Hamiltonian, bipartite graph with radius 3,
        diameter 3, and girth 4. ::

            sage: G.is_hamiltonian()
            True
            sage: G.is_bipartite()
            True
            sage: G.radius()
            3
            sage: G.diameter()
            3
            sage: G.girth()
            4

        It is a perfect, triangle-free graph having chromatic number 2. ::

            sage: G.is_perfect()
            True
            sage: G.is_triangle_free()
            True
            sage: G.chromatic_number()
            2
        """
        edge_dict = {
            0: [1,5,6],
            1: [2,7],
            2: [3,8],
            3: [4,9],
            4: [5,10],
            5: [11],
            6: [7,9],
            7: [10],
            8: [9,11],
            10: [11]}
        pos_dict = {
            0: [2, 0],
            1: [1, 1.73205080756888],
            2: [-1, 1.73205080756888],
            3: [-2, 0],
            4: [-1, -1.73205080756888],
            5: [1, -1.73205080756888],
            6: [1, 0],
            7: [0.5, 0.866025403784439],
            8: [-0.5, 0.866025403784439],
            9: [-1, 0],
            10: [-0.5, -0.866025403784439],
            11: [0.5, -0.866025403784439]}
        return graph.Graph(edge_dict, pos=pos_dict, name="Franklin graph")

    def FruchtGraph(self):
        """
        Returns a Frucht Graph.

        A Frucht graph has 12 nodes and 18 edges. It is the smallest cubic
        identity graph. It is planar and it is Hamiltonian.

        This constructor is dependent on NetworkX's numeric labeling.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the first
        seven nodes are on the outer circle, with the next four on an inner
        circle and the last in the center.

        REFERENCES:

        - [1] Weisstein, E. (1999). "Frucht Graph - from Wolfram
          MathWorld". [Online] Available:
          http://mathworld.wolfram.com/FruchtGraph.html [2007, February 17]

        EXAMPLES::

            sage: FRUCHT = graphs.FruchtGraph()
            sage: FRUCHT
            Frucht graph: Graph on 12 vertices
            sage: FRUCHT.graph6_string()
            'KhCKM?_EGK?L'
            sage: (graphs.FruchtGraph()).show() # long time
        """
        pos_dict = {}
        for i in range(7):
            x = float(2*(cos((pi/2) + ((2*pi)/7)*i)))
            y = float(2*(sin((pi/2) + ((2*pi)/7)*i)))
            pos_dict[i] = (x,y)
        pos_dict[7] = (0,1)
        pos_dict[8] = (-1,0)
        pos_dict[9] = (0,-1)
        pos_dict[10] = (1,0)
        pos_dict[11] = (0,0)
        import networkx
        G = networkx.frucht_graph()
        return graph.Graph(G, pos=pos_dict, name="Frucht graph")

    def GoldnerHararyGraph(self):
        r"""
        Return the Goldner-Harary graph.

        For more information, see this
        `Wikipedia article on the Goldner-Harary graph <http://en.wikipedia.org/wiki/Goldner%E2%80%93Harary_graph>`_.

        EXAMPLES:

        The Goldner-Harary graph is named after A. Goldner and Frank Harary.
        It is a planar graph having 11 vertices and 27 edges. ::

            sage: G = graphs.GoldnerHararyGraph(); G
            Goldner-Harary graph: Graph on 11 vertices
            sage: G.is_planar()
            True
            sage: G.order()
            11
            sage: G.size()
            27

        The Goldner-Harary graph is chordal with radius 2, diameter 2, and
        girth 3. ::

            sage: G.is_chordal()
            True
            sage: G.radius()
            2
            sage: G.diameter()
            2
            sage: G.girth()
            3

        Its chromatic number is 4 and its automorphism group is isomorphic to
        the dihedral group `D_6`. ::

            sage: G.chromatic_number()
            4
            sage: ag = G.automorphism_group()
            sage: ag.is_isomorphic(DihedralGroup(6))
            True
        """
        edge_dict = {
            0: [1,3,4],
            1: [2,3,4,5,6,7,10],
            2: [3,7],
            3: [7,8,9,10],
            4: [3,5,9,10],
            5: [10],
            6: [7,10],
            7: [8,10],
            8: [10],
            9: [10]}

        pos = {
            0: (-2, 0),
            1: (0, 1.5),
            2: (2, 0),
            3: (0, -1.5),
            4: (-1.5, 0),
            5: (-0.5, 0.5),
            6: (0.5, 0.5),
            7: (1.5, 0),
            8: (0.5, -0.5),
            9: (-0.5, -0.5),
            10: (0, 0)}

        return graph.Graph(edge_dict, pos = pos, name="Goldner-Harary graph")

    def GrayGraph(self, embedding=1):
        r"""
        Returns the Gray graph.

        See the :wikipedia:`Wikipedia page on the Gray Graph
        <Gray_graph>`.

        INPUT:

        - ``embedding`` -- two embeddings are available, and can be
          selected by setting ``embedding`` to 1 or 2.

        EXAMPLES::

            sage: g = graphs.GrayGraph()
            sage: g.order()
            54
            sage: g.size()
            81
            sage: g.girth()
            8
            sage: g.diameter()
            6
            sage: g.show(figsize=[10, 10])   # long time
            sage: graphs.GrayGraph(embedding = 2).show(figsize=[10, 10])   # long time

        TESTS::

            sage: graphs.GrayGraph(embedding = 3)
            Traceback (most recent call last):
            ...
            ValueError: The value of embedding must be 1, 2, or 3.
        """

        g = graphs.LCFGraph(54, [-25,7,-7,13,-13,25], 9)
        g.name("Gray graph")

        if embedding == 1:
            o = g.automorphism_group(orbits=True)[-1]
            _circle_embedding(g, o[0], center=(0, 0), radius=1)
            _circle_embedding(g, o[1], center=(0, 0), radius=.6, shift=-.5)

        elif embedding != 2:
            raise ValueError("The value of embedding must be 1, 2, or 3.")

        return g

    def GrotzschGraph(self):
        r"""
        Returns the Grötzsch graph.

        The Grötzsch graph is an example of a triangle-free graph with
        chromatic number equal to 4. For more information, see this
        `Wikipedia article on Grötzsch graph <http://en.wikipedia.org/wiki/Gr%C3%B6tzsch_graph>`_.

        REFERENCE:

        - [1] Weisstein, Eric W. "Grotzsch Graph."
          From MathWorld--A Wolfram Web Resource.
          http://mathworld.wolfram.com/GroetzschGraph.html

        EXAMPLES:

        The Grötzsch graph is named after Herbert Grötzsch. It is a
        Hamiltonian graph with 11 vertices and 20 edges. ::

            sage: G = graphs.GrotzschGraph(); G
            Grotzsch graph: Graph on 11 vertices
            sage: G.is_hamiltonian()
            True
            sage: G.order()
            11
            sage: G.size()
            20

        The Grötzsch graph is triangle-free and having radius 2, diameter 2,
        and girth 4. ::

            sage: G.is_triangle_free()
            True
            sage: G.radius()
            2
            sage: G.diameter()
            2
            sage: G.girth()
            4

        Its chromatic number is 4 and its automorphism group is isomorphic
        to the dihedral group `D_5`. ::

            sage: G.chromatic_number()
            4
            sage: ag = G.automorphism_group()
            sage: ag.is_isomorphic(DihedralGroup(5))
            True
        """
        g = graph.Graph()
        g.add_vertices(range(11))

        edges = [];
        for u in range(1,6):
            edges.append( (0,u) )

        edges.append( (10,6) )

        for u in range(6,10):
            edges.append( (u,u+1) )
            edges.append( (u,u-4) )

        edges.append( (10,1) )

        for u in range(7,11):
            edges.append( (u,u-6) )

        edges.append((6,5))

        g.add_edges(edges)

        pos = {}
        pos[0] = (0,0)
        for u in range(1,6):
            theta = (u-1)*2*pi/5
            pos[u] = (float(5*sin(theta)),float(5*cos(theta)))
            pos[u+5] = (2*pos[u][0], 2*pos[u][1])

        g.set_pos(pos)
        g.name("Grotzsch graph")
        return g

    def HeawoodGraph(self):
        """
        Returns a Heawood graph.

        The Heawood graph is a cage graph that has 14 nodes. It is a cubic
        symmetric graph. (See also the Moebius-Kantor graph). It is
        nonplanar and Hamiltonian. It has diameter = 3, radius = 3, girth =
        6, chromatic number = 2. It is 4-transitive but not 5-transitive.

        This constructor is dependent on NetworkX's numeric labeling.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the nodes are
        positioned in a circular layout with the first node appearing at
        the top, and then continuing counterclockwise.

        REFERENCES:

        - [1] Weisstein, E. (1999). "Heawood Graph - from Wolfram
          MathWorld". [Online] Available:
          http://mathworld.wolfram.com/HeawoodGraph.html [2007, February 17]

        EXAMPLES::

            sage: H = graphs.HeawoodGraph()
            sage: H
            Heawood graph: Graph on 14 vertices
            sage: H.graph6_string()
            'MhEGHC@AI?_PC@_G_'
            sage: (graphs.HeawoodGraph()).show() # long time
        """
        pos_dict = {}
        for i in range(14):
            x = float(cos((pi/2) + (pi/7)*i))
            y = float(sin((pi/2) + (pi/7)*i))
            pos_dict[i] = (x,y)
        import networkx
        G = networkx.heawood_graph()
        return graph.Graph(G, pos=pos_dict, name="Heawood graph")

    def HerschelGraph(self):
        r"""
        Returns the Herschel graph.

        For more information, see this
        `Wikipedia article on the Herschel graph <http://en.wikipedia.org/wiki/Herschel_graph>`_.

        EXAMPLES:

        The Herschel graph is named after Alexander Stewart Herschel. It is
        a planar, bipartite graph with 11 vertices and 18 edges. ::

            sage: G = graphs.HerschelGraph(); G
            Herschel graph: Graph on 11 vertices
            sage: G.is_planar()
            True
            sage: G.is_bipartite()
            True
            sage: G.order()
            11
            sage: G.size()
            18

        The Herschel graph is a perfect graph with radius 3, diameter 4, and
        girth 4. ::

            sage: G.is_perfect()
            True
            sage: G.radius()
            3
            sage: G.diameter()
            4
            sage: G.girth()
            4

        Its chromatic number is 2 and its automorphism group is
        isomorphic to the dihedral group `D_6`. ::

            sage: G.chromatic_number()
            2
            sage: ag = G.automorphism_group()
            sage: ag.is_isomorphic(DihedralGroup(6))
            True
        """
        edge_dict = {
            0: [1,3,4],
            1: [2,5,6],
            2: [3,7],
            3: [8,9],
            4: [5,9],
            5: [10],
            6: [7,10],
            7: [8],
            8: [10],
            9: [10]}
        pos_dict = {
            0: [2, 0],
            1: [0, 2],
            2: [-2, 0],
            3: [0, -2],
            4: [1, 0],
            5: [0.5, 0.866025403784439],
            6: [-0.5, 0.866025403784439],
            7: [-1, 0],
            8: [-0.5, -0.866025403784439],
            9: [0.5, -0.866025403784439],
            10: [0, 0]}
        return graph.Graph(edge_dict, pos=pos_dict, name="Herschel graph")

    def HigmanSimsGraph(self, relabel=True):
        r"""
        The Higman-Sims graph is a remarkable strongly regular
        graph of degree 22 on 100 vertices.  For example, it can
        be split into two sets of 50 vertices each, so that each
        half induces a subgraph isomorphic to the
        Hoffman-Singleton graph
        (:meth:`~HoffmanSingletonGraph`).
        This can be done in 352 ways (see [BROUWER-HS-2009]_).

        Its most famous property is that the automorphism
        group has an index 2 subgroup which is one of the
        26 sporadic groups. [HIGMAN1968]_

        The construction used here follows [HAFNER2004]_.

        INPUT:

        - ``relabel`` - default: ``True``.  If ``True`` the
          vertices will be labeled with consecutive integers.
          If ``False`` the labels are strings that are three
          digits long. "xyz" means the vertex is in group
          x (zero through three), pentagon or pentagram y
          (zero through four), and is vertex z (zero
          through four) of that pentagon or pentagram.
          See [HAFNER2004]_ for more.

        OUTPUT:

        The Higman-Sims graph.

        EXAMPLES:

        A split into the first 50 and last 50 vertices
        will induce two copies of the Hoffman-Singleton graph,
        and we illustrate another such split, which is obvious
        based on the construction used. ::

            sage: H = graphs.HigmanSimsGraph()
            sage: A = H.subgraph(range(0,50))
            sage: B = H.subgraph(range(50,100))
            sage: K = graphs.HoffmanSingletonGraph()
            sage: K.is_isomorphic(A) and K.is_isomorphic(B)
            True
            sage: C = H.subgraph(range(25,75))
            sage: D = H.subgraph(range(0,25)+range(75,100))
            sage: K.is_isomorphic(C) and K.is_isomorphic(D)
            True

        The automorphism group contains only one nontrivial
        proper normal subgroup, which is of index 2 and is
        simple.  It is known as the Higman-Sims group.  ::

            sage: H = graphs.HigmanSimsGraph()
            sage: G = H.automorphism_group()
            sage: g=G.order(); g
            88704000
            sage: K = G.normal_subgroups()[1]
            sage: K.is_simple()
            True
            sage: g//K.order()
            2

        REFERENCES:

            .. [BROUWER-HS-2009] `Higman-Sims graph
               <http://www.win.tue.nl/~aeb/graphs/Higman-Sims.html>`_.
               Andries E. Brouwer, accessed 24 October 2009.
            .. [HIGMAN1968] A simple group of order 44,352,000,
               Math.Z. 105 (1968) 110-113. D.G. Higman & C. Sims.
            .. [HAFNER2004] `On the graphs of Hoffman-Singleton and
               Higman-Sims
               <http://www.combinatorics.org/Volume_11/PDF/v11i1r77.pdf>`_.
               The Electronic Journal of Combinatorics 11 (2004), #R77,
               Paul R. Hafner, accessed 24 October 2009.

        AUTHOR:

            - Rob Beezer (2009-10-24)
        """
        HS = graph.Graph()
        HS.name('Higman-Sims graph')

        # Four groups of either five pentagons, or five pentagrams
        # 4 x 5 x 5 = 100 vertices
        # First digit is "group", second is "penta{gon|gram}", third is "vertex"
        vlist = ['%d%d%d'%(g,p,v)
                        for g in range(4) for p in range(5) for v in range(5)]
        for avertex in vlist:
            HS.add_vertex(avertex)

        # Edges: Within groups 0 and 2, joined as pentagons
        # Edges: Within groups 1 and 3, joined as pentagrams
        for g in range(4):
            shift = 1
            if g in [1,3]:
                shift += 1
            for p in range(5):
                for v in range(5):
                    HS.add_edge(('%d%d%d'%(g,p,v), '%d%d%d'%(g,p,(v+shift)%5)))

        # Edges: group 0 to group 1
        for x in range(5):
            for m in range(5):
                for c in range(5):
                    y = (m*x+c)%5
                    HS.add_edge(('0%d%d'%(x,y), '1%d%d'%(m,c)))

        # Edges: group 1 to group 2
        for m in range(5):
            for A in range(5):
                for B in range(5):
                    c = (2*(m-A)*(m-A)+B)%5
                    HS.add_edge(('1%d%d'%(m,c), '2%d%d'%(A,B)))

        # Edges: group 2 to group 3
        for A in range(5):
            for a in range(5):
                for b in range(5):
                    B = (2*A*A+3*a*A-a*a+b)%5
                    HS.add_edge(('2%d%d'%(A,B), '3%d%d'%(a,b)))

        # Edges: group 3 to group 0
        for a in range(5):
            for b in range(5):
                for x in range(5):
                    y = ((x-a)*(x-a)+b)%5
                    HS.add_edge(('3%d%d'%(a,b), '0%d%d'%(x,y)))

        # Edges: group 0 to group 2
        for x in range(5):
            for A in range(5):
                for B in range(5):
                    y = (3*x*x+A*x+B+1)%5
                    HS.add_edge(('0%d%d'%(x,y), '2%d%d'%(A,B)))
                    y = (3*x*x+A*x+B-1)%5
                    HS.add_edge(('0%d%d'%(x,y), '2%d%d'%(A,B)))

        # Edges: group 1 to group 3
        for m in range(5):
            for a in range(5):
                for b in range(5):
                    c = (m*(m-a)+b+2)%5
                    HS.add_edge(('1%d%d'%(m,c), '3%d%d'%(a,b)))
                    c = (m*(m-a)+b-2)%5
                    HS.add_edge(('1%d%d'%(m,c), '3%d%d'%(a,b)))

        # Rename to integer vertex labels, creating dictionary
        # Or not, and create identity mapping
        if relabel:
            vmap = HS.relabel(return_map=True)
        else:
            vmap={}
            for v in vlist:
                vmap[v] = v
        # Layout vertices in a circle
        # In the order given in vlist
        # Using labels from vmap
        pos_dict = {}
        for i in range(100):
            x = float(cos((pi/2) + ((2*pi)/100)*i))
            y = float(sin((pi/2) + ((2*pi)/100)*i))
            pos_dict[vmap[vlist[i]]] = (x,y)
        HS.set_pos(pos_dict)
        return HS

    def HoffmanSingletonGraph(self):
        r"""
        Returns the Hoffman-Singleton graph.

        The Hoffman-Singleton graph is the Moore graph of degree 7,
        diameter 2 and girth 5. The Hoffman-Singleton theorem states that
        any Moore graph with girth 5 must have degree 2, 3, 7 or 57. The
        first three respectively are the pentagon, the Petersen graph, and
        the Hoffman-Singleton graph. The existence of a Moore graph with
        girth 5 and degree 57 is still open.

        A Moore graph is a graph with diameter `d` and girth
        `2d + 1`. This implies that the graph is regular, and
        distance regular.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. A novel algorithm written by
        Tom Boothby gives a random layout which is pleasing to the eye.

        REFERENCES:

        - [1] Godsil, C. and Royle, G. Algebraic Graph Theory.
          Springer, 2001.

        EXAMPLES::

            sage: HS = graphs.HoffmanSingletonGraph()
            sage: Set(HS.degree())
            {7}
            sage: HS.girth()
            5
            sage: HS.diameter()
            2
            sage: HS.num_verts()
            50

        Note that you get a different layout each time you create the graph.
        ::

            sage: HS.layout()[1]
            (-0.844..., 0.535...)
            sage: graphs.HoffmanSingletonGraph().layout()[1]
            (-0.904..., 0.425...)

        """
        H = graph.Graph({ \
        'q00':['q01'], 'q01':['q02'], 'q02':['q03'], 'q03':['q04'], 'q04':['q00'], \
        'q10':['q11'], 'q11':['q12'], 'q12':['q13'], 'q13':['q14'], 'q14':['q10'], \
        'q20':['q21'], 'q21':['q22'], 'q22':['q23'], 'q23':['q24'], 'q24':['q20'], \
        'q30':['q31'], 'q31':['q32'], 'q32':['q33'], 'q33':['q34'], 'q34':['q30'], \
        'q40':['q41'], 'q41':['q42'], 'q42':['q43'], 'q43':['q44'], 'q44':['q40'], \
        'p00':['p02'], 'p02':['p04'], 'p04':['p01'], 'p01':['p03'], 'p03':['p00'], \
        'p10':['p12'], 'p12':['p14'], 'p14':['p11'], 'p11':['p13'], 'p13':['p10'], \
        'p20':['p22'], 'p22':['p24'], 'p24':['p21'], 'p21':['p23'], 'p23':['p20'], \
        'p30':['p32'], 'p32':['p34'], 'p34':['p31'], 'p31':['p33'], 'p33':['p30'], \
        'p40':['p42'], 'p42':['p44'], 'p44':['p41'], 'p41':['p43'], 'p43':['p40']})
        for j in range(5):
            for i in range(5):
                for k in range(5):
                    con = (i+j*k)%5
                    H.add_edge(('q%d%d'%(k,con),'p%d%d'%(j,i)))
        H.name('Hoffman-Singleton graph')
        from sage.combinat.combinat import permutations
        from sage.misc.prandom import randint
        P = permutations([1,2,3,4])
        qpp = [0]+P[randint(0,23)]
        ppp = [0]+P[randint(0,23)]
        qcycle = lambda i,s : ['q%s%s'%(i,(j+s)%5) for j in qpp]
        pcycle = lambda i,s : ['p%s%s'%(i,(j+s)%5) for j in ppp]
        l = 0
        s = 0
        D = []
        while l < 5:
            for q in qcycle(l,s):
                D.append(q)
            vv = 'p%s'%q[1]
            s = int([v[-1] for v in H.neighbors(q) if v[:2] == vv][0])
            for p in pcycle(l,s):
                D.append(p)
            vv = 'q%s'%(int(p[1])+1)
            v = [v[-1] for v in H.neighbors(p) if v[:2] == vv]
            if len(v):
                s = int(v[0])
            l+=1
        map = H.relabel(return_map=True)
        pos_dict = {}
        for i in range(50):
            x = float(cos((pi/2) + ((2*pi)/50)*i))
            y = float(sin((pi/2) + ((2*pi)/50)*i))
            pos_dict[map[D[i]]] = (x,y)
        H.set_pos(pos_dict)
        return H

    def HoffmanGraph(self):
        r"""
        Returns the Hoffman Graph.

        See the :wikipedia:`Wikipedia page on the Hoffman graph
        <Hoffman_graph>`.

        EXAMPLES::

            sage: g = graphs.HoffmanGraph()
            sage: g.is_bipartite()
            True
            sage: g.is_hamiltonian() # long time
            True
            sage: g.radius()
            3
            sage: g.diameter()
            4
            sage: g.automorphism_group().cardinality()
            48
        """
        g = graph.Graph({
                0: [1, 7, 8, 13],
                1: [2, 9, 14],
                2: [3, 8, 10],
                3: [4, 9, 15],
                4: [5, 10, 11],
                5: [6, 12, 14],
                6: [7, 11, 13],
                7: [12, 15],
                8: [12, 14],
                9: [11, 13],
                10: [12, 15],
                11: [14],
                13: [15]})
        g.set_pos({})
        _circle_embedding(g, range(8))
        _circle_embedding(g, range(8, 14), radius=.7, shift=.5)
        _circle_embedding(g, [14, 15], radius=.1)

        g.name("Hoffman Graph")

        return g

    def HoltGraph(self):
        r"""
        Returns the Holt Graph.

        See the :wikipedia:`Wikipedia page on the Holt graph
        <Holt_graph>`.

        EXAMPLES::

            sage: g = graphs.HoltGraph()
            sage: g.chromatic_number()
            3
            sage: g.is_hamiltonian() # long time
            True
            sage: g.radius()
            3
            sage: g.diameter()
            4
            sage: g.girth()
            5
            sage: g.automorphism_group().cardinality()
            18
        """
        g = graph.Graph({
                0: [9, 12],
                1: [11, 14],
                2: [13, 16],
                3: [15, 18],
                4: [17, 20],
                5: [19, 22],
                6: [21, 24],
                7: [23, 26],
                8: [10, 25]
                },pos={})

        g.add_vertices(range(27))
        g.add_cycle(range(9))

        g.add_cycle([13,21,11,19,9,17,25,15,23])
        g.add_cycle([12,16,20, 24, 10, 14, 18, 22, 26])

        _circle_embedding(g, range(9), shift = .75)
        _circle_embedding(g, range(9, 27), radius = .7, shift = 0)

        g.name("Holt graph")

        return g

    def LjubljanaGraph(self, embedding=1):
        r"""
        Returns the Ljubljana Graph.

        The Ljubljana graph is a bipartite 3-regular graph on 112
        vertices and 168 edges. It is not vertex-transitive as it has
        two orbits which are also independent sets of size 56. See the
        :wikipedia:`Wikipedia page on the Ljubljana Graph
        <Ljubljana_graph>`.

        The default embedding is obtained from the Heawood graph.

        INPUT:

        - ``embedding`` -- two embeddings are available, and can be
          selected by setting ``embedding`` to 1 or 2.

        EXAMPLES::

            sage: g = graphs.LjubljanaGraph()
            sage: g.order()
            112
            sage: g.size()
            168
            sage: g.girth()
            10
            sage: g.diameter()
            8
            sage: g.show(figsize=[10, 10])   # long time
            sage: graphs.LjubljanaGraph(embedding=2).show(figsize=[10, 10])   # long time

        TESTS::

            sage: graphs.LjubljanaGraph(embedding=3)
            Traceback (most recent call last):
            ...
            ValueError: The value of embedding must be 1 or 2.
        """

        L = [47, -23, -31, 39, 25, -21, -31, -41, 25, 15, 29, -41, -19, 15,
             -49, 33, 39, -35, -21, 17, -33, 49, 41, 31, -15, -29, 41, 31,
             -15, -25, 21, 31, -51, -25, 23, 9, -17, 51, 35, -29, 21, -51,
             -39, 33, -9, -51, 51, -47, -33, 19, 51, -21, 29, 21, -31, -39]

        g = graphs.LCFGraph(112, L, 2)
        g.name("Ljubljana graph")

        if embedding == 1:
            return g

        elif embedding == 2:
            dh = graphs.HeawoodGraph().get_pos()

            # Correspondence between the vertices of the Heawood Graph and
            # 8-sets of the Ljubljana Graph.

            d = {
                0: [1, 21, 39, 57, 51, 77, 95, 107],
                1: [2, 22, 38, 58, 50, 78, 94, 106],
                2: [3, 23, 37, 59, 49, 79, 93, 105],
                3: [4, 24, 36, 60, 48, 80, 92, 104],
                4: [5, 25, 35, 61, 15, 81, 91, 71],
                9: [6, 26, 44, 62, 16, 82, 100, 72],
                10: [7, 27, 45, 63, 17, 83, 101, 73],
                11: [8, 28, 46, 64, 18, 84, 102, 74],
                12: [9, 29, 47, 65, 19, 85, 103, 75],
                13: [10, 30, 0, 66, 20, 86, 56, 76],
                8: [11, 31, 111, 67, 99, 87, 55, 43],
                7: [12, 32, 110, 68, 98, 88, 54, 42],
                6: [13, 33, 109, 69, 97, 89, 53, 41],
                5: [14, 34, 108, 70, 96, 90, 52, 40]
                }

            # The vertices of each 8-set are plotted on a circle, and the
            # circles are slowly shifted to obtain a symmetric drawing.

            for i, (u, vertices) in enumerate(d.iteritems()):
                _circle_embedding(g, vertices, center=dh[u], radius=.1,
                        shift=8.*i/14)

            return g

        else:
            raise ValueError("The value of embedding must be 1 or 2.")

    def KneserGraph(self,n,k):
        r"""
        Returns the Kneser Graph with parameters `n, k`.

        The Kneser Graph with parameters `n,k` is the graph
        whose vertices are the `k`-subsets of `[0,1,\dots,n-1]`, and such
        that two vertices are adjacent if their corresponding sets
        are disjoint.

        For example, the Petersen Graph can be defined
        as the Kneser Graph with parameters `5,2`.

        EXAMPLE::

            sage: KG=graphs.KneserGraph(5,2)
            sage: print KG.vertices()
            [{4, 5}, {1, 3}, {2, 5}, {2, 3}, {3, 4}, {3, 5}, {1, 4}, {1, 5}, {1, 2}, {2, 4}]
            sage: P=graphs.PetersenGraph()
            sage: P.is_isomorphic(KG)
            True

        TESTS::

            sage: KG=graphs.KneserGraph(0,0)
            Traceback (most recent call last):
            ...
            ValueError: Parameter n should be a strictly positive integer
            sage: KG=graphs.KneserGraph(5,6)
            Traceback (most recent call last):
            ...
            ValueError: Parameter k should be a strictly positive integer inferior to n
        """

        if not n>0:
            raise ValueError, "Parameter n should be a strictly positive integer"
        if not (k>0 and k<=n):
            raise ValueError, "Parameter k should be a strictly positive integer inferior to n"

        g=graph.Graph(name="Kneser graph with parameters "+str(n)+","+str(k))
        from sage.combinat.subset import Subsets

        if k>n/2:
            g.add_vertices(Subsets(n,k).list())

        S = Subsets(n,k)
        for s in S:
            for t in Subsets(S.s.difference(s),k):
                g.add_edge(s,t)

        return g

    def McGeeGraph(self, embedding=2):
        r"""
        Returns the McGee Graph.

        See the :wikipedia:`Wikipedia page on the McGee Graph
        <McGee_graph>`.

        INPUT:

        - ``embedding`` -- two embeddings are available, and can be
          selected by setting ``embedding`` to 1 or 2.

        EXAMPLES::

            sage: g = graphs.McGeeGraph()
            sage: g.order()
            24
            sage: g.size()
            36
            sage: g.girth()
            7
            sage: g.diameter()
            4
            sage: g.show()
            sage: graphs.McGeeGraph(embedding=1).show()

        TESTS::

            sage: graphs.McGeeGraph(embedding=3)
            Traceback (most recent call last):
            ...
            ValueError: The value of embedding must be 1 or 2.
        """

        L = [47, -23, -31, 39, 25, -21, -31, -41, 25, 15, 29, -41, -19, 15,
             -49, 33, 39, -35, -21, 17, -33, 49, 41, 31, -15, -29, 41, 31,
             -15, -25, 21, 31, -51, -25, 23, 9, -17, 51, 35, -29, 21, -51,
             -39, 33, -9, -51, 51, -47, -33, 19, 51, -21, 29, 21, -31, -39]

        g = graphs.LCFGraph(24, [12, 7, -7], 8)
        g.name('McGee graph')

        if embedding == 1:
            return g

        elif embedding == 2:

            o = [[7, 2, 13, 8, 19, 14, 1, 20],
                 [5, 4, 11, 10, 17, 16, 23, 22],
                 [3, 12, 9, 18, 15, 0, 21, 6]]

            _circle_embedding(g, o[0], radius=1.5)
            _circle_embedding(g, o[1], radius=3, shift=-.5)
            _circle_embedding(g, o[2], radius=2.25, shift=.5)

            return g

        else:
            raise ValueError("The value of embedding must be 1 or 2.")

    def MoebiusKantorGraph(self):
        """
        Returns a Moebius-Kantor Graph.

        A Moebius-Kantor graph is a cubic symmetric graph. (See also the
        Heawood graph). It has 16 nodes and 24 edges. It is nonplanar and
        Hamiltonian. It has diameter = 4, girth = 6, and chromatic number =
        2. It is identical to the Generalized Petersen graph, P[8,3].

        PLOTTING: See the plotting section for the generalized Petersen graphs.

        REFERENCES:

        - [1] Weisstein, E. (1999). "Moebius-Kantor Graph - from
          Wolfram MathWorld". [Online] Available:
          http://mathworld.wolfram.com/Moebius-KantorGraph.html [2007,
          February 17]

        EXAMPLES::

            sage: MK = graphs.MoebiusKantorGraph()
            sage: MK
            Moebius-Kantor Graph: Graph on 16 vertices
            sage: MK.graph6_string()
            'OhCGKE?O@?ACAC@I?Q_AS'
            sage: (graphs.MoebiusKantorGraph()).show() # long time
        """
        G=graphs.GeneralizedPetersenGraph(8,3)
        G.name("Moebius-Kantor Graph")
        return G

    def MoserSpindle(self):
        r"""
        Returns the Moser spindle.

        For more information, see this
        `MathWorld article on the Moser spindle <http://mathworld.wolfram.com/MoserSpindle.html>`_.

        EXAMPLES:

        The Moser spindle is a planar graph having 7 vertices and 11 edges. ::

            sage: G = graphs.MoserSpindle(); G
            Moser spindle: Graph on 7 vertices
            sage: G.is_planar()
            True
            sage: G.order()
            7
            sage: G.size()
            11

        It is a Hamiltonian graph with radius 2, diameter 2, and girth 3. ::

            sage: G.is_hamiltonian()
            True
            sage: G.radius()
            2
            sage: G.diameter()
            2
            sage: G.girth()
            3

        The Moser spindle has chromatic number 4 and its automorphism
        group is isomorphic to the dihedral group `D_4`. ::

            sage: G.chromatic_number()
            4
            sage: ag = G.automorphism_group()
            sage: ag.is_isomorphic(DihedralGroup(4))
            True
        """
        edge_dict = {
            0: [1,4,5,6],
            1: [2,5],
            2: [3,5],
            3: [4,6],
            4: [6]}
        pos_dict = {
            0: [0, 2],
            1: [-1.90211303259031, 0.618033988749895],
            2: [-1.17557050458495, -1.61803398874989],
            3: [1.17557050458495, -1.61803398874989],
            4: [1.90211303259031, 0.618033988749895],
            5: [1, 0],
            6: [-1, 0]}
        return graph.Graph(edge_dict, pos=pos_dict, name="Moser spindle")

    def MycielskiGraph(self, k=1, relabel=True):
        r"""
        Returns the `k`-th Mycielski Graph.

        The graph `M_k` is triangle-free and has chromatic number
        equal to `k`. These graphs show, constructively, that there
        are triangle-free graphs with arbitrarily high chromatic
        number.

        The Mycielski graphs are built recursively starting with
        `M_0`, an empty graph; `M_1`, a single vertex graph; and `M_2`
        is the graph `K_2`.  `M_{k+1}` is then built from `M_k`
        as follows:

        If the vertices of `M_k` are `v_1,\ldots,v_n`, then the
        vertices of `M_{k+1}` are
        `v_1,\ldots,v_n,w_1,\ldots,w_n,z`. Vertices `v_1,\ldots,v_n`
        induce a copy of `M_k`. Vertices `w_1,\ldots,w_n` are an
        independent set. Vertex `z` is adjacent to all the
        `w_i`-vertices. Finally, vertex `w_i` is adjacent to vertex
        `v_j` iff `v_i` is adjacent to `v_j`.

        INPUT:

        - ``k`` Number of steps in the construction process.

        - ``relabel`` Relabel the vertices so their names are the integers
          ``range(n)`` where ``n`` is the number of vertices in the graph.

        EXAMPLE:

        The Mycielski graph `M_k` is triangle-free and has chromatic
        number equal to `k`. ::

            sage: g = graphs.MycielskiGraph(5)
            sage: g.is_triangle_free()
            True
            sage: g.chromatic_number()
            5

        The graphs `M_4` is (isomorphic to) the Grotzsch graph. ::

            sage: g = graphs.MycielskiGraph(4)
            sage: g.is_isomorphic(graphs.GrotzschGraph())
            True

        REFERENCES:

        -  [1] Weisstein, Eric W. "Mycielski Graph."
           From MathWorld--A Wolfram Web Resource.
           http://mathworld.wolfram.com/MycielskiGraph.html

        """
        g = graph.Graph()
        g.name("Mycielski Graph " + str(k))

        if k<0:
            raise ValueError, "parameter k must be a nonnegative integer"

        if k == 0:
            return g

        if k == 1:
            g.add_vertex(0)
            return g

        if k == 2:
            g.add_edge(0,1)
            return g

        g0 = graphs.MycielskiGraph(k-1)
        g = graphs.MycielskiStep(g0)
        g.name("Mycielski Graph " + str(k))
        if relabel: g.relabel()

        return g

    def MycielskiStep(self, g):
        r"""
        Perform one iteration of the Mycielski construction.

        See the documentation for ``MycielskiGraph`` which uses this
        method. We expose it to all users in case they may find it
        useful.

        EXAMPLE. One iteration of the Mycielski step applied to the
        5-cycle yields a graph isomorphic to the Grotzsch graph ::

            sage: g = graphs.CycleGraph(5)
            sage: h = graphs.MycielskiStep(g)
            sage: h.is_isomorphic(graphs.GrotzschGraph())
            True
        """

        # Make a copy of the input graph g
        gg = g.copy()

        # rename a vertex v of gg as (1,v)
        renamer = dict( [ (v, (1,v)) for v in g.vertices() ] )
        gg.relabel(renamer)

        # add the w vertices to gg as (2,v)
        wlist = [ (2,v) for v in g.vertices() ]
        gg.add_vertices(wlist)

        # add the z vertex as (0,0)
        gg.add_vertex((0,0))

        # add the edges from z to w_i
        gg.add_edges( [ ( (0,0) , (2,v) ) for v in g.vertices() ] )

        # make the v_i w_j edges
        for v in g.vertices():
            gg.add_edges( [ ((1,v),(2,vv)) for vv in g.neighbors(v) ] )

        return gg

    def NauruGraph(self, embedding=2):
        """
        Returns the Nauru Graph.

        See the :wikipedia:`Wikipedia page on the Nauru Graph
        <Nauru_graph>`.

        INPUT:

        - ``embedding`` -- two embeddings are available, and can be
          selected by setting ``embedding`` to 1 or 2.

        EXAMPLES::

            sage: g = graphs.NauruGraph()
            sage: g.order()
            24
            sage: g.size()
            36
            sage: g.girth()
            6
            sage: g.diameter()
            4
            sage: g.show()
            sage: graphs.NauruGraph(embedding=1).show()

        TESTS::

            sage: graphs.NauruGraph(embedding=3)
            Traceback (most recent call last):
            ...
            ValueError: The value of embedding must be 1 or 2.
            sage: graphs.NauruGraph(embedding=1).is_isomorphic(g)
            True
        """

        if embedding == 1:
            g = graphs.LCFGraph(24, [5, -9, 7, -7, 9, -5], 4)
            g.name('Nauru Graph')
            return g
        elif embedding == 2:
            g = graphs.GeneralizedPetersenGraph(12, 5)
            g.name("Nauru Graph")
            return g
        else:
            raise ValueError("The value of embedding must be 1 or 2.")

    def PappusGraph(self):
        """
        Returns the Pappus graph, a graph on 18 vertices.

        The Pappus graph is cubic, symmetric, and distance-regular.

        EXAMPLES::

            sage: G = graphs.PappusGraph()
            sage: G.show()  # long time
            sage: L = graphs.LCFGraph(18, [5,7,-7,7,-7,-5], 3)
            sage: L.show()  # long time
            sage: G.is_isomorphic(L)
            True
        """
        pos_dict = {}
        for i in range(6):
            pos_dict[i] = [float(cos(pi/2 + ((2*pi)/6)*i)),\
                           float(sin(pi/2 + ((2*pi)/6)*i))]
            pos_dict[6 + i] = [(2/3.0)*float(cos(pi/2 + ((2*pi)/6)*i)),\
                               (2/3.0)*float(sin(pi/2 + ((2*pi)/6)*i))]
            pos_dict[12 + i] = [(1/3.0)*float(cos(pi/2 + ((2*pi)/6)*i)),\
                                (1/3.0)*float(sin(pi/2 + ((2*pi)/6)*i))]
        return graph.Graph({0:[1,5,6],1:[2,7],2:[3,8],3:[4,9],4:[5,10],\
                            5:[11],6:[13,17],7:[12,14],8:[13,15],9:[14,16],\
                            10:[15,17],11:[12,16],12:[15],13:[16],14:[17]},\
                           pos=pos_dict, name="Pappus Graph")

    def PetersenGraph(self):
        """
        The Petersen Graph is a named graph that consists of 10 vertices
        and 15 edges, usually drawn as a five-point star embedded in a
        pentagon.

        The Petersen Graph is a common counterexample. For example, it is
        not Hamiltonian.

        PLOTTING: See the plotting section for the generalized Petersen graphs.

        EXAMPLES: We compare below the Petersen graph with the default
        spring-layout versus a planned position dictionary of [x,y]
        tuples::

            sage: petersen_spring = Graph({0:[1,4,5], 1:[0,2,6], 2:[1,3,7], 3:[2,4,8], 4:[0,3,9], 5:[0,7,8], 6:[1,8,9], 7:[2,5,9], 8:[3,5,6], 9:[4,6,7]})
            sage: petersen_spring.show() # long time
            sage: petersen_database = graphs.PetersenGraph()
            sage: petersen_database.show() # long time
        """
        P=graphs.GeneralizedPetersenGraph(5,2)
        P.name("Petersen graph")
        return P

    def ShrikhandeGraph(self):
        """
        Returns the Shrikhande graph.

        For more information, see the `MathWorld article on the Shrikhande graph
        <http://mathworld.wolfram.com/ShrikhandeGraph.html>`_ or the `Wikipedia
        article on the Shrikhande graph
        <http://en.wikipedia.org/wiki/Shrikhande_graph>`_.

        EXAMPLES:

        The Shrikhande graph was defined by S. S. Shrikhande in 1959. It has
        `16` vertices and `48` edges, and is strongly regular of degree `6` with
        parameters `(2,2)`::

            sage: G = graphs.ShrikhandeGraph(); G
            Shrikhande graph: Graph on 16 vertices
            sage: G.order()
            16
            sage: G.size()
            48
            sage: G.is_regular(6)
            True
            sage: set([ len([x for x in G.neighbors(i) if x in G.neighbors(j)])
            ...     for i in range(G.order())
            ...     for j in range(i) ])
            set([2])

        It is non-planar, and both Hamiltonian and Eulerian::

            sage: G.is_planar()
            False
            sage: G.is_hamiltonian()
            True
            sage: G.is_eulerian()
            True

        It has radius `2`, diameter `2`, and girth `3`::

            sage: G.radius()
            2
            sage: G.diameter()
            2
            sage: G.girth()
            3

        Its chromatic number is `4` and its automorphism group is of order
        `192`::

            sage: G.chromatic_number()
            4
            sage: G.automorphism_group().cardinality()
            192

        It is an integral graph since it has only integral eigenvalues::

            sage: G.characteristic_polynomial().factor()
            (x - 6) * (x - 2)^6 * (x + 2)^9

        It is a toroidal graph, and its embedding on a torus is dual to an
        embedding of the Dyck graph (:meth:`DyckGraph <GraphGenerators.DyckGraph>`).
        """
        pos_dict = {}
        for i in range(8):
            pos_dict[i] = [float(cos((2*i) * pi/8)),
                           float(sin((2*i) * pi/8))]
            pos_dict[8 + i] = [0.5 * pos_dict[i][0],
                               0.5 * pos_dict[i][1]]
        edge_dict = {
            0O00: [0O06, 0O07, 0O01, 0O02,   0O11, 0O17],
            0O01: [0O07, 0O00, 0O02, 0O03,   0O12, 0O10],
            0O02: [0O00, 0O01, 0O03, 0O04,   0O13, 0O11],
            0O03: [0O01, 0O02, 0O04, 0O05,   0O14, 0O12],
            0O04: [0O02, 0O03, 0O05, 0O06,   0O15, 0O13],
            0O05: [0O03, 0O04, 0O06, 0O07,   0O16, 0O14],
            0O06: [0O04, 0O05, 0O07, 0O00,   0O17, 0O15],
            0O07: [0O05, 0O06, 0O00, 0O01,   0O10, 0O16],

            0O10: [0O12, 0O13, 0O15, 0O16,   0O07, 0O01],
            0O11: [0O13, 0O14, 0O16, 0O17,   0O00, 0O02],
            0O12: [0O14, 0O15, 0O17, 0O10,   0O01, 0O03],
            0O13: [0O15, 0O16, 0O10, 0O11,   0O02, 0O04],
            0O14: [0O16, 0O17, 0O11, 0O12,   0O03, 0O05],
            0O15: [0O17, 0O10, 0O12, 0O13,   0O04, 0O06],
            0O16: [0O10, 0O11, 0O13, 0O14,   0O05, 0O07],
            0O17: [0O11, 0O12, 0O14, 0O15,   0O06, 0O00]
        }

        return graph.Graph(edge_dict, pos=pos_dict, name="Shrikhande graph")

    def ThomsenGraph(self):
        """
        Returns the Thomsen Graph.

        The Thomsen Graph is actually a complete bipartite graph with (n1,
        n2) = (3, 3). It is also called the Utility graph.

        PLOTTING: See CompleteBipartiteGraph.

        EXAMPLES::

            sage: T = graphs.ThomsenGraph()
            sage: T
            Thomsen graph: Graph on 6 vertices
            sage: T.graph6_string()
            'EFz_'
            sage: (graphs.ThomsenGraph()).show() # long time
        """
        pos_dict = {0:(-1,1),1:(0,1),2:(1,1),3:(-1,0),4:(0,0),5:(1,0)}
        import networkx
        G = networkx.complete_bipartite_graph(3,3)
        return graph.Graph(G, pos=pos_dict, name="Thomsen graph")

    def Tutte12Cage(self):
        r"""
        Returns Tutte's 12-Cage.

        See the :wikipedia:`Wikipedia page on the Tutte 12-Cage
        <Tutte_12-cage>`.

        EXAMPLES::

            sage: g = graphs.Tutte12Cage()
            sage: g.order()
            126
            sage: g.size()
            189
            sage: g.girth()
            12
            sage: g.diameter()
            6
            sage: g.show()
        """
        L = [17, 27, -13, -59, -35, 35, -11, 13, -53, 53, -27, 21, 57, 11,
             -21, -57, 59, -17]

        g = graphs.LCFGraph(126, L, 7)
        g.name("Tutte 12-Cage")
        return g

    def TutteCoxeterGraph(self, embedding=2):
        r"""
        Returns the Tutte-Coxeter graph.

        See the :wikipedia:`Wikipedia page on the Tutte-Coxeter Graph
        <Tutte-Coxeter_graph>`.

        INPUT:

        - ``embedding`` -- two embeddings are available, and can be
          selected by setting ``embedding`` to 1 or 2.

        EXAMPLES::

            sage: g = graphs.TutteCoxeterGraph()
            sage: g.order()
            30
            sage: g.size()
            45
            sage: g.girth()
            8
            sage: g.diameter()
            4
            sage: g.show()
            sage: graphs.TutteCoxeterGraph(embedding=1).show()

        TESTS::

            sage: graphs.TutteCoxeterGraph(embedding=3)
            Traceback (most recent call last):
            ...
            ValueError: The value of embedding must be 1 or 2.
        """

        g = graphs.LCFGraph(30, [-13, -9, 7, -7, 9, 13], 5)
        g.name("Tutte-Coxeter graph")

        if embedding == 1:
            d = {
                0: [1, 3, 5, 7, 29],
                1: [2, 4, 6, 28, 0],
                2: [8, 18, 26, 22, 12],
                3: [9, 13, 23, 27, 17],
                4: [11, 15, 21, 25, 19],
                5: [10, 14, 24, 20, 16]
                }

            _circle_embedding(g, d[0], center=(-1, 1), radius=.25)
            _circle_embedding(g, d[1], center=(1, 1), radius=.25)
            _circle_embedding(g, d[2], center=(-.8, 0), radius=.25, shift=2.5)
            _circle_embedding(g, d[3], center=(1.2, 0), radius=.25)
            _circle_embedding(g, d[4], center=(-1, -1), radius=.25, shift=2)
            _circle_embedding(g, d[5], center=(1, -1), radius=.25)

            return g

        elif embedding == 2:
            return g

        else:
            raise ValueError("The value of embedding must be 1 or 2.")

    def WagnerGraph(self):
        """
        Returns the Wagner Graph.

        See the :wikipedia:`Wikipedia page on the Wagner Graph
        <Wagner_graph>`.

        EXAMPLES::

            sage: g = graphs.WagnerGraph()
            sage: g.order()
            8
            sage: g.size()
            12
            sage: g.girth()
            4
            sage: g.diameter()
            2
            sage: g.show()
        """
        g = graphs.LCFGraph(8, [4], 8)
        g.name("Wagner Graph")
        return g

###########################################################################
#   Families of Graphs
###########################################################################

    def BalancedTree(self, r, h):
        r"""
        Returns the perfectly balanced tree of height `h \geq 1`,
        whose root has degree `r \geq 2`.

        The number of vertices of this graph is
        `1 + r + r^2 + \cdots + r^h`, that is,
        `\frac{r^{h+1} - 1}{r - 1}`. The number of edges is one
        less than the number of vertices.

        INPUT:

        - ``r`` -- positive integer `\geq 2`. The degree of the root node.

        - ``h`` -- positive integer `\geq 1`. The height of the balanced tree.

        OUTPUT:

        The perfectly balanced tree of height `h \geq 1` and whose root has
        degree `r \geq 2`. A ``NetworkXError`` is returned if `r < 2` or
        `h < 1`.

        ALGORITHM:

        Uses `NetworkX <http://networkx.lanl.gov>`_.

        EXAMPLES:

        A balanced tree whose root node has degree `r = 2`, and of height
        `h = 1`, has order 3 and size 2::

            sage: G = graphs.BalancedTree(2, 1); G
            Balanced tree: Graph on 3 vertices
            sage: G.order(); G.size()
            3
            2
            sage: r = 2; h = 1
            sage: v = 1 + r
            sage: v; v - 1
            3
            2

        Plot a balanced tree of height 5, whose root node has degree `r = 3`::

            sage: G = graphs.BalancedTree(3, 5)
            sage: G.show()   # long time

        A tree is bipartite. If its vertex set is finite, then it is planar. ::

            sage: r = randint(2, 5); h = randint(1, 7)
            sage: T = graphs.BalancedTree(r, h)
            sage: T.is_bipartite()
            True
            sage: T.is_planar()
            True
            sage: v = (r^(h + 1) - 1) / (r - 1)
            sage: T.order() == v
            True
            sage: T.size() == v - 1
            True

        TESTS:

         Normally we would only consider balanced trees whose root node
         has degree `r \geq 2`, but the construction degenerates
         gracefully::

            sage: graphs.BalancedTree(1, 10)
            Balanced tree: Graph on 2 vertices

            sage: graphs.BalancedTree(-1, 10)
            Balanced tree: Graph on 1 vertex

        Similarly, we usually want the tree must have height `h \geq 1`
        but the algorithm also degenerates gracefully here::

            sage: graphs.BalancedTree(3, 0)
            Balanced tree: Graph on 1 vertex

            sage: graphs.BalancedTree(5, -2)
            Balanced tree: Graph on 0 vertices

            sage: graphs.BalancedTree(-2,-2)
            Balanced tree: Graph on 0 vertices
        """
        import networkx
        return graph.Graph(networkx.balanced_tree(r, h), name="Balanced tree")

    def BubbleSortGraph(self, n):
        r"""
        Returns the bubble sort graph `B(n)`.

        The vertices of the bubble sort graph are the set of permutations on
        `n` symbols. Two vertices are adjacent if one can be obtained from the
        other by swapping the labels in the `i`-th and `(i+1)`-th positions for
        `1 \leq i \leq n-1`. In total, `B(n)` has order `n!`. Thus, the order
        of `B(n)` increases according to `f(n) = n!`.

        INPUT:

        - ``n`` -- positive integer. The number of symbols to permute.

        OUTPUT:

        The bubble sort graph `B(n)` on `n` symbols. If `n < 1`, a
        ``ValueError`` is returned.

        EXAMPLES::

            sage: g = graphs.BubbleSortGraph(4); g
            Bubble sort: Graph on 24 vertices
            sage: g.plot() # long time

        The bubble sort graph on `n = 1` symbol is the trivial graph `K_1`::

            sage: graphs.BubbleSortGraph(1)
            Bubble sort: Graph on 1 vertex

        If `n \geq 1`, then the order of `B(n)` is `n!`::

            sage: n = randint(1, 8)
            sage: g = graphs.BubbleSortGraph(n)
            sage: g.order() == factorial(n)
            True

        TESTS:

        Input ``n`` must be positive::

            sage: graphs.BubbleSortGraph(0)
            Traceback (most recent call last):
            ...
            ValueError: Invalid number of symbols to permute, n should be >= 1
            sage: graphs.BubbleSortGraph(randint(-10^6, 0))
            Traceback (most recent call last):
            ...
            ValueError: Invalid number of symbols to permute, n should be >= 1

        AUTHORS:

        - Michael Yurko (2009-09-01)
        """
        # sanity checks
        if n < 1:
            raise ValueError(
                "Invalid number of symbols to permute, n should be >= 1")
        if n == 1:
            return graph.Graph(self.CompleteGraph(n), name="Bubble sort")
        from sage.combinat.permutation import Permutations
        #create set from which to permute
        label_set = [str(i) for i in xrange(1, n + 1)]
        d = {}
        #iterate through all vertices
        for v in Permutations(label_set):
            tmp_dict = {}
            #add all adjacencies
            for i in xrange(n - 1):
                #swap entries
                v[i], v[i + 1] = v[i + 1], v[i]
                #add new vertex
                new_vert = ''.join(v)
                tmp_dict[new_vert] = None
                #swap back
                v[i], v[i + 1] = v[i + 1], v[i]
            #add adjacency dict
            d[''.join(v)] = tmp_dict
        return graph.Graph(d, name="Bubble sort")

    def CirculantGraph(self, n, adjacency):
        r"""
        Returns a circulant graph with n nodes.

        A circulant graph has the property that the vertex i is connected
        with the vertices i+j and i-j for each j in adj.

        INPUT:


        -  ``n`` - number of vertices in the graph

        -  ``adjacency`` - the list of j values


        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, each circulant
        graph will be displayed with the first (0) node at the top, with
        the rest following in a counterclockwise manner.

        Filling the position dictionary in advance adds O(n) to the
        constructor.

        EXAMPLES: Compare plotting using the predefined layout and
        networkx::

            sage: import networkx
            sage: n = networkx.cycle_graph(23)
            sage: spring23 = Graph(n)
            sage: posdict23 = graphs.CirculantGraph(23,2)
            sage: spring23.show() # long time
            sage: posdict23.show() # long time

        We next view many cycle graphs as a Sage graphics array. First we
        use the ``CirculantGraph`` constructor, which fills in
        the position dictionary::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.CirculantGraph(i+3,i)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time

        Compare to plotting with the spring-layout algorithm::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    spr = networkx.cycle_graph(i+3)
            ...    k = Graph(spr)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time

        Passing a 1 into adjacency should give the cycle.

        ::

            sage: graphs.CirculantGraph(6,1)==graphs.CycleGraph(6)
            True
            sage: graphs.CirculantGraph(7,[1,3]).edges(labels=false)
            [(0, 1),
            (0, 3),
            (0, 4),
            (0, 6),
            (1, 2),
            (1, 4),
            (1, 5),
            (2, 3),
            (2, 5),
            (2, 6),
            (3, 4),
            (3, 6),
            (4, 5),
            (5, 6)]
        """
        if not isinstance(adjacency,list):
            adjacency=[adjacency]
        pos_dict = {}
        for i in range(n):
            x = float(cos((pi/2) + ((2*pi)/n)*i))
            y = float(sin((pi/2) + ((2*pi)/n)*i))
            pos_dict[i] = (x,y)
        G=graph.Graph(n, name="Circulant graph ("+str(adjacency)+")")
        G._pos=pos_dict
        for v in G:
            G.add_edges([(v,(v+j)%n) for j in adjacency])
            G.add_edges([(v,(v-j)%n) for j in adjacency])
        return G

    def CompleteGraph(self, n):
        """
        Returns a complete graph on n nodes.

        A Complete Graph is a graph in which all nodes are connected to all
        other nodes.

        This constructor is dependent on vertices numbered 0 through n-1 in
        NetworkX complete_graph()

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, each complete
        graph will be displayed with the first (0) node at the top, with
        the rest following in a counterclockwise manner.

        In the complete graph, there is a big difference visually in using
        the spring-layout algorithm vs. the position dictionary used in
        this constructor. The position dictionary flattens the graph,
        making it clear which nodes an edge is connected to. But the
        complete graph offers a good example of how the spring-layout
        works. The edges push outward (everything is connected), causing
        the graph to appear as a 3-dimensional pointy ball. (See examples
        below).

        EXAMPLES: We view many Complete graphs with a Sage Graphics Array,
        first with this constructor (i.e., the position dictionary
        filled)::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.CompleteGraph(i+3)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time

        We compare to plotting with the spring-layout algorithm::

            sage: import networkx
            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    spr = networkx.complete_graph(i+3)
            ...    k = Graph(spr)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time

        Compare the constructors (results will vary)

        ::

            sage: import networkx
            sage: t = cputime()
            sage: n = networkx.complete_graph(389); spring389 = Graph(n)
            sage: cputime(t)           # random
            0.59203700000000126
            sage: t = cputime()
            sage: posdict389 = graphs.CompleteGraph(389)
            sage: cputime(t)           # random
            0.6680419999999998

        We compare plotting::

            sage: import networkx
            sage: n = networkx.complete_graph(23)
            sage: spring23 = Graph(n)
            sage: posdict23 = graphs.CompleteGraph(23)
            sage: spring23.show() # long time
            sage: posdict23.show() # long time
        """
        pos_dict = {}
        for i in range(n):
            x = float(cos((pi/2) + ((2*pi)/n)*i))
            y = float(sin((pi/2) + ((2*pi)/n)*i))
            pos_dict[i] = (x,y)
        import networkx
        G = networkx.complete_graph(n)
        return graph.Graph(G, pos=pos_dict, name="Complete graph")

    def CompleteBipartiteGraph(self, n1, n2):
        """
        Returns a Complete Bipartite Graph sized n1+n2, with each of the
        nodes [0,(n1-1)] connected to each of the nodes [n1,(n2-1)] and
        vice versa.

        A Complete Bipartite Graph is a graph with its vertices partitioned
        into two groups, V1 and V2. Each v in V1 is connected to every v in
        V2, and vice versa.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, each complete
        bipartite graph will be displayed with the first n1 nodes on the
        top row (at y=1) from left to right. The remaining n2 nodes appear
        at y=0, also from left to right. The shorter row (partition with
        fewer nodes) is stretched to the same length as the longer row,
        unless the shorter row has 1 node; in which case it is centered.
        The x values in the plot are in domain [0,maxn1,n2].

        In the Complete Bipartite graph, there is a visual difference in
        using the spring-layout algorithm vs. the position dictionary used
        in this constructor. The position dictionary flattens the graph and
        separates the partitioned nodes, making it clear which nodes an
        edge is connected to. The Complete Bipartite graph plotted with the
        spring-layout algorithm tends to center the nodes in n1 (see
        spring_med in examples below), thus overlapping its nodes and
        edges, making it typically hard to decipher.

        Filling the position dictionary in advance adds O(n) to the
        constructor. Feel free to race the constructors below in the
        examples section. The much larger difference is the time added by
        the spring-layout algorithm when plotting. (Also shown in the
        example below). The spring model is typically described as
        `O(n^3)`, as appears to be the case in the NetworkX source
        code.

        EXAMPLES: Two ways of constructing the complete bipartite graph,
        using different layout algorithms::

            sage: import networkx
            sage: n = networkx.complete_bipartite_graph(389,157); spring_big = Graph(n)   # long time
            sage: posdict_big = graphs.CompleteBipartiteGraph(389,157)                    # long time

        Compare the plotting::

            sage: n = networkx.complete_bipartite_graph(11,17)
            sage: spring_med = Graph(n)
            sage: posdict_med = graphs.CompleteBipartiteGraph(11,17)

        Notice here how the spring-layout tends to center the nodes of n1

        ::

            sage: spring_med.show() # long time
            sage: posdict_med.show() # long time

        View many complete bipartite graphs with a Sage Graphics Array,
        with this constructor (i.e., the position dictionary filled)::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.CompleteBipartiteGraph(i+1,4)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time

        We compare to plotting with the spring-layout algorithm::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    spr = networkx.complete_bipartite_graph(i+1,4)
            ...    k = Graph(spr)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time

        Trac ticket #12155::

            sage: graphs.CompleteBipartiteGraph(5,6).complement()
            complement(Complete bipartite graph): Graph on 11 vertices
        """
        pos_dict = {}
        c1 = 1 # scaling factor for top row
        c2 = 1 # scaling factor for bottom row
        c3 = 0 # pad to center if top row has 1 node
        c4 = 0 # pad to center if bottom row has 1 node
        if n1 > n2:
            if n2 == 1:
                c4 = (n1-1)/2
            else:
                c2 = ((n1-1)/(n2-1))
        elif n2 > n1:
            if n1 == 1:
                c3 = (n2-1)/2
            else:
                c1 = ((n2-1)/(n1-1))
        for i in range(n1):
            x = c1*i + c3
            y = 1
            pos_dict[i] = (x,y)
        for i in range(n1+n2)[n1:]:
            x = c2*(i-n1) + c4
            y = 0
            pos_dict[i] = (x,y)
        import networkx
        from sage.graphs.graph import Graph
        G = networkx.complete_bipartite_graph(n1,n2)
        return Graph(G, pos=pos_dict, name="Complete bipartite graph")

    def CompleteMultipartiteGraph(self, l):
        r"""
        Returns a complete multipartite graph.

        INPUT:

        - ``l`` -- a list of integers : the respective sizes
          of the components.

        EXAMPLE:

        A complete tripartite graph with sets of sizes
        `5, 6, 8`::

            sage: g = graphs.CompleteMultipartiteGraph([5, 6, 8]); g
            Multipartite Graph with set sizes [5, 6, 8]: Graph on 19 vertices

        It clearly has a chromatic number of 3::

            sage: g.chromatic_number()
            3
        """

        from sage.graphs.graph import Graph
        g = Graph()
        for i in l:
            g = g + self.CompleteGraph(i)

        g = g.complement()
        g.name("Multipartite Graph with set sizes "+str(l))

        return g

    def CubeGraph(self, n):
        r"""
        Returns the hypercube in `n` dimensions.

        The hypercube in `n` dimension is build upon the binary
        strings on `n` bits, two of them being adjacent if
        they differ in exactly one bit. Hence, the distance
        between two vertices in the hypercube is the Hamming
        distance.

        EXAMPLES:

        The distance between `0100110` and `1011010` is
        `5`, as expected ::

            sage: g = graphs.CubeGraph(7)
            sage: g.distance('0100110','1011010')
            5

        Plot several `n`-cubes in a Sage Graphics Array ::

            sage: g = []
            sage: j = []
            sage: for i in range(6):
            ...    k = graphs.CubeGraph(i+1)
            ...    g.append(k)
            ...
            sage: for i in range(2):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show(figsize=[6,4]) # long time

        Use the plot options to display larger `n`-cubes

        ::

            sage: g = graphs.CubeGraph(9)
            sage: g.show(figsize=[12,12],vertex_labels=False, vertex_size=20) # long time

        AUTHORS:

        - Robert Miller
        """
        theta = float(pi/n)

        d = {'':[]}
        dn={}
        p = {'':(float(0),float(0))}
        pn={}

        # construct recursively the adjacency dict and the positions
        for i in xrange(n):
            ci = float(cos(i*theta))
            si = float(sin(i*theta))
            for v,e in d.iteritems():
                v0 = v+'0'
                v1 = v+'1'
                l0 = [v1]
                l1 = [v0]
                for m in e:
                    l0.append(m+'0')
                    l1.append(m+'1')
                dn[v0] = l0
                dn[v1] = l1
                x,y = p[v]
                pn[v0] = (x, y)
                pn[v1] = (x+ci, y+si)
            d,dn = dn,{}
            p,pn = pn,{}

        # construct the graph
        r = graph.Graph(name="%d-Cube"%n)
        r.add_vertices(d.keys())
        for u,L in d.iteritems():
            for v in L:
                r.add_edge(u,v)
        r.set_pos(p)

        return r

    def FriendshipGraph(self, n):
        r"""
        Returns the friendship graph `F_n`.

        The friendship graph is also known as the Dutch windmill graph. Let
        `C_3` be the cycle graph on 3 vertices. Then `F_n` is constructed by
        joining `n \geq 1` copies of `C_3` at a common vertex. If `n = 1`,
        then `F_1` is isomorphic to `C_3` (the triangle graph). If `n = 2`,
        then `F_2` is the butterfly graph, otherwise known as the bowtie
        graph. For more information, see this
        `Wikipedia article on the friendship graph <http://en.wikipedia.org/wiki/Friendship_graph>`_.

        INPUT:

        - ``n`` -- positive integer; the number of copies of `C_3` to use in
          constructing `F_n`.

        OUTPUT:

        - The friendship graph `F_n` obtained from `n` copies of the cycle
          graph `C_3`.

        .. seealso::

            - :meth:`GraphGenerators.ButterflyGraph`

        EXAMPLES:

        The first few friendship graphs. ::

            sage: A = []; B = []
            sage: for i in range(9):
            ...       g = graphs.FriendshipGraph(i + 1)
            ...       A.append(g)
            sage: for i in range(3):
            ...       n = []
            ...       for j in range(3):
            ...           n.append(A[3*i + j].plot(vertex_size=20, vertex_labels=False))
            ...       B.append(n)
            sage: G = sage.plot.graphics.GraphicsArray(B)
            sage: G.show()  # long time

        For `n = 1`, the friendship graph `F_1` is isomorphic to the cycle
        graph `C_3`, whose visual representation is a triangle. ::

            sage: G = graphs.FriendshipGraph(1); G
            Friendship graph: Graph on 3 vertices
            sage: G.show()  # long time
            sage: G.is_isomorphic(graphs.CycleGraph(3))
            True

        For `n = 2`, the friendship graph `F_2` is isomorphic to the
        butterfly graph, otherwise known as the bowtie graph. ::

            sage: G = graphs.FriendshipGraph(2); G
            Friendship graph: Graph on 5 vertices
            sage: G.is_isomorphic(graphs.ButterflyGraph())
            True

        If `n \geq 1`, then the friendship graph `F_n` has `2n + 1` vertices
        and `3n` edges. It has radius 1, diameter 2, girth 3, and
        chromatic number 3. Furthermore, `F_n` is planar and Eulerian. ::

            sage: n = randint(1, 10^3)
            sage: G = graphs.FriendshipGraph(n)
            sage: G.order() == 2*n + 1
            True
            sage: G.size() == 3*n
            True
            sage: G.radius()
            1
            sage: G.diameter()
            2
            sage: G.girth()
            3
            sage: G.chromatic_number()
            3
            sage: G.is_planar()
            True
            sage: G.is_eulerian()
            True

        TESTS:

        The input ``n`` must be a positive integer. ::

            sage: graphs.FriendshipGraph(randint(-10^5, 0))
            Traceback (most recent call last):
            ...
            ValueError: n must be a positive integer
        """
        # sanity checks
        if n < 1:
            raise ValueError("n must be a positive integer")
        # construct the friendship graph
        if n == 1:
            G = self.CycleGraph(3)
            G.name("Friendship graph")
            return G
        # build the edge and position dictionaries
        from sage.functions.trig import cos, sin
        from sage.rings.real_mpfr import RR
        from sage.symbolic.constants import pi
        N = 2*n + 1           # order of F_n
        d = (2*pi) / (N - 1)  # angle between external nodes
        edge_dict = {}
        pos_dict = {}
        for i in range(N - 2):
            if i & 1:  # odd numbered node
                edge_dict.setdefault(i, [i + 1, N - 1])
            else:      # even numbered node
                edge_dict.setdefault(i, [N - 1])
            pos_dict.setdefault(i, [RR(cos(i*d)), RR(sin(i*d))])
        edge_dict.setdefault(N - 2, [0, N - 1])
        pos_dict.setdefault(N - 2, [RR(cos(d * (N-2))), RR(sin(d * (N-2)))])
        pos_dict.setdefault(N - 1, [0, 0])
        return graph.Graph(edge_dict, pos=pos_dict, name="Friendship graph")

    def FuzzyBallGraph(self, partition, q):
        r"""
        Construct a Fuzzy Ball graph with the integer partition
        ``partition`` and ``q`` extra vertices.

        Let `q` be an integer and let `m_1,m_2,...,m_k` be a set of positive
        integers.  Let `n=q+m_1+...+m_k`.  The Fuzzy Ball graph with partition
        `m_1,m_2,...,m_k` and `q` extra vertices is the graph constructed from the
        graph `G=K_n` by attaching, for each `i=1,2,...,k`, a new vertex `a_i` to
        `m_i` distinct vertices of `G`.

        For given positive integers `k` and `m` and nonnegative
        integer `q`, the set of graphs ``FuzzyBallGraph(p, q)`` for
        all partitions `p` of `m` with `k` parts are cospectral with
        respect to the normalized Laplacian.

        EXAMPLES::

            sage: graphs.FuzzyBallGraph([3,1],2).adjacency_matrix()
            [0 1 1 1 1 1 1 0]
            [1 0 1 1 1 1 1 0]
            [1 1 0 1 1 1 1 0]
            [1 1 1 0 1 1 0 1]
            [1 1 1 1 0 1 0 0]
            [1 1 1 1 1 0 0 0]
            [1 1 1 0 0 0 0 0]
            [0 0 0 1 0 0 0 0]


        Pick positive integers `m` and `k` and a nonnegative integer `q`.
        All the FuzzyBallGraphs constructed from partitions of `m` with
        `k` parts should be cospectral with respect to the normalized
        Laplacian::

            sage: m=4; q=2; k=2
            sage: g_list=[graphs.FuzzyBallGraph(p,q) for p in Partitions(m, length=k)]
            sage: set([g.laplacian_matrix(normalized=True).charpoly() for g in g_list])  # long time (7s on sage.math, 2011)
            set([x^8 - 8*x^7 + 4079/150*x^6 - 68689/1350*x^5 + 610783/10800*x^4 - 120877/3240*x^3 + 1351/100*x^2 - 931/450*x])
        """
        if len(partition)<1:
            raise ValueError, "partition must be a nonempty list of positive integers"
        n=q+sum(partition)
        g=graphs.CompleteGraph(n)
        curr_vertex=0
        for e,p in enumerate(partition):
            g.add_edges([(curr_vertex+i, 'a{0}'.format(e+1)) for i in range(p)])
            curr_vertex+=p
        return g

    def FibonacciTree(self, n):
        r"""
        Returns the graph of the Fibonacci Tree `F_{i}` of order `n`.
        `F_{i}` is recursively defined as the a tree with a root vertex
        and two attached child trees `F_{i-1}` and `F_{i-2}`, where
        `F_{1}` is just one vertex and `F_{0}` is empty.

        INPUT:

        - ``n`` - the recursion depth of the Fibonacci Tree

        EXAMPLES::

            sage: g = graphs.FibonacciTree(3)
            sage: g.is_tree()
            True

        ::

            sage: l1 = [ len(graphs.FibonacciTree(_)) + 1 for _ in range(6) ]
            sage: l2 = list(fibonacci_sequence(2,8))
            sage: l1 == l2
            True

        AUTHORS:

        - Harald Schilly and Yann Laigle-Chapuy (2010-03-25)
        """
        T = graph.Graph(name="Fibonacci-Tree-%d"%n)
        if n == 1: T.add_vertex(0)
        if n < 2: return T

        from sage.combinat.combinat import fibonacci_sequence
        F = list(fibonacci_sequence(n + 2))
        s = 1.618 ** (n / 1.618 - 1.618)
        pos = {}

        def fib(level, node, y):
            pos[node] = (node, y)
            if level < 2: return
            level -= 1
            y -= s
            diff = F[level]
            T.add_edge(node, node - diff)
            if level == 1: # only one child
                pos[node - diff] = (node, y)
                return
            T.add_edge(node, node + diff)
            fib(level, node - diff, y)
            fib(level - 1, node + diff, y)

        T.add_vertices(xrange(sum(F[:-1])))
        fib(n, F[n + 1] - 1, 0)
        T.set_pos(pos)

        return T

    def GeneralizedPetersenGraph(self, n,k):
        r"""
        Returns a generalized Petersen graph with `2n` nodes. The variables
        `n`, `k` are integers such that `n>2` and `0<k\leq\lfloor(n-1)`/`2\rfloor`

        For `k=1` the result is a graph isomorphic to the circular ladder graph
        with the same `n`. The regular Petersen Graph has `n=5` and `k=2`.
        Other named graphs that can be described using this notation include
        the Desargues graph and the Moebius-Kantor graph.

        INPUT:

        - ``n`` - the number of nodes is `2*n`.

        - ``k`` - integer `0<k\leq\lfloor(n-1)`/`2\rfloor`. Decides
          how inner vertices are connected.

        PLOTTING: Upon construction, the position dictionary is filled to
        override the spring-layout algorithm. By convention, the generalized
        Petersen graphs are displayed as an inner and outer cycle pair, with
        the first n nodes drawn on the outer circle. The first (0) node is
        drawn at the top of the outer-circle, moving counterclockwise after that.
        The inner circle is drawn with the (n)th node at the top, then
        counterclockwise as well.

        EXAMPLES: For `k=1` the resulting graph will be isomorphic to a circular
        ladder graph. ::

            sage: g = graphs.GeneralizedPetersenGraph(13,1)
            sage: g2 = graphs.CircularLadderGraph(13)
            sage: g.is_isomorphic(g2)
            True

        The Desargues graph::

            sage: g = graphs.GeneralizedPetersenGraph(10,3)
            sage: g.girth()
            6
            sage: g.is_bipartite()
            True

        AUTHORS:

        - Anders Jonsson (2009-10-15)
        """
        if (n < 3):
                raise ValueError("n must be larger than 2")
        if (k < 1 or k>((n-1)/2)):
                raise ValueError("k must be in 1<= k <=floor((n-1)/2)")
        pos_dict = {}
        G=graphs.EmptyGraph()
        for i in range(n):
            x = float(cos((pi/2) + ((2*pi)/n)*i))
            y = float(sin((pi/2) + ((2*pi)/n)*i))
            pos_dict[i] = (x,y)
        for i in range(n, 2*n):
            x = float(0.5*cos((pi/2) + ((2*pi)/n)*i))
            y = float(0.5*sin((pi/2) + ((2*pi)/n)*i))
            pos_dict[i] = (x,y)
        for i in range(n):
            G.add_edge(i, (i+1) % n)
            G.add_edge(i, i+n)
            G.add_edge(i+n, n + (i+k) % n)
        return graph.Graph(G, pos=pos_dict, name="Generalized Petersen graph (n="+str(n)+",k="+str(k)+")")

    def HyperStarGraph(self,n,k):
        r"""
        Returns the hyper-star graph HS(n,k).

        The vertices of the hyper-star graph are the set of binary strings
        of length n which contain k 1s. Two vertices, u and v, are adjacent
        only if u can be obtained from v by swapping the first bit with a
        different symbol in another position.

        INPUT:

        -  ``n``

        -  ``k``

        EXAMPLES::

            sage: g = graphs.HyperStarGraph(6,3)
            sage: g.plot() # long time

        REFERENCES:

        - Lee, Hyeong-Ok, Jong-Seok Kim, Eunseuk Oh, and Hyeong-Seok Lim.
          "Hyper-Star Graph: A New Interconnection Network Improving the
          Network Cost of the Hypercube." In Proceedings of the First EurAsian
          Conference on Information and Communication Technology, 858-865.
          Springer-Verlag, 2002.

        AUTHORS:

        - Michael Yurko (2009-09-01)
        """
        from sage.combinat.combination import Combinations
        # dictionary associating the positions of the 1s to the corresponding
        # string: e.g. if n=6 and k=3, comb_to_str([0,1,4])=='110010'
        comb_to_str={}
        for c in Combinations(n,k):
            L = ['0']*n
            for i in c:
                L[i]='1'
            comb_to_str[tuple(c)] = ''.join(L)

        g=graph.Graph(name="HS(%d,%d)"%(n,k))
        g.add_vertices(comb_to_str.values())

        for c in Combinations(range(1,n),k): # 0 is not in c
            L = []
            u = comb_to_str[tuple(c)]
            # switch 0 with the 1s
            for i in xrange(len(c)):
                v = tuple([0]+c[:i]+c[i+1:])
                g.add_edge( u , comb_to_str[v] )

        return g

    def LCFGraph(self, n, shift_list, repeats):
        """
        Returns the cubic graph specified in LCF notation.

        LCF (Lederberg-Coxeter-Fruchte) notation is a concise way of
        describing cubic Hamiltonian graphs. The way a graph is constructed
        is as follows. Since there is a Hamiltonian cycle, we first create
        a cycle on n nodes. The variable shift_list = [s_0, s_1, ...,
        s_k-1] describes edges to be created by the following scheme: for
        each i, connect vertex i to vertex (i + s_i). Then, repeats
        specifies the number of times to repeat this process, where on the
        jth repeat we connect vertex (i + j\*len(shift_list)) to vertex (
        i + j\*len(shift_list) + s_i).

        INPUT:


        -  ``n`` - the number of nodes.

        -  ``shift_list`` - a list of integer shifts mod n.

        -  ``repeats`` - the number of times to repeat the
           process.


        EXAMPLES::

            sage: G = graphs.LCFGraph(4, [2,-2], 2)
            sage: G.is_isomorphic(graphs.TetrahedralGraph())
            True

        ::

            sage: G = graphs.LCFGraph(20, [10,7,4,-4,-7,10,-4,7,-7,4], 2)
            sage: G.is_isomorphic(graphs.DodecahedralGraph())
            True

        ::

            sage: G = graphs.LCFGraph(14, [5,-5], 7)
            sage: G.is_isomorphic(graphs.HeawoodGraph())
            True

        The largest cubic nonplanar graph of diameter three::

            sage: G = graphs.LCFGraph(20, [-10,-7,-5,4,7,-10,-7,-4,5,7,-10,-7,6,-5,7,-10,-7,5,-6,7], 1)
            sage: G.degree()
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
            sage: G.diameter()
            3
            sage: G.show()  # long time

        PLOTTING: LCF Graphs are plotted as an n-cycle with edges in the
        middle, as described above.

        REFERENCES:

        - [1] Frucht, R. "A Canonical Representation of Trivalent
          Hamiltonian Graphs." J. Graph Th. 1, 45-60, 1976.

        - [2] Grunbaum, B.  Convex Polytope es. New York: Wiley,
          pp. 362-364, 1967.

        - [3] Lederberg, J. 'DENDRAL-64: A System for Computer
          Construction, Enumeration and Notation of Organic Molecules
          as Tree Structures and Cyclic Graphs. Part II. Topology of
          Cyclic Graphs.' Interim Report to the National Aeronautics
          and Space Administration. Grant NsG 81-60. December 15,
          1965.  http://profiles.nlm.nih.gov/BB/A/B/I/U/_/bbabiu.pdf.
        """
        import networkx
        pos_dict = {}
        for i in range(n):
            x = float(cos(pi/2 + ((2*pi)/n)*i))
            y = float(sin(pi/2 + ((2*pi)/n)*i))
            pos_dict[i] = [x,y]
        return graph.Graph(networkx.LCF_graph(n, shift_list, repeats),\
                           pos=pos_dict, name="LCF Graph")

    def NKStarGraph(self,n,k):
        r"""
        Returns the (n,k)-star graph.

        The vertices of the (n,k)-star graph are the set of all arrangements of
        n symbols into labels of length k. There are two adjacency rules for
        the (n,k)-star graph. First, two vertices are adjacent if one can be
        obtained from the other by swapping the first symbol with another
        symbol. Second, two vertices are adjacent if one can be obtained from
        the other by swapping the first symbol with an external symbol (a
        symbol not used in the original label).

        INPUT:

        -  ``n``

        -  ``k``

        EXAMPLES::

            sage: g = graphs.NKStarGraph(4,2)
            sage: g.plot() # long time

        REFERENCES:

        - Wei-Kuo, Chiang, and Chen Rong-Jaye. "The (n, k)-star graph: A
          generalized star graph." Information Processing Letters 56,
          no. 5 (December 8, 1995): 259-264.

        AUTHORS:

        - Michael Yurko (2009-09-01)
        """
        from sage.combinat.permutation import Arrangements
        #set from which to permute
        set = [str(i) for i in xrange(1,n+1)]
        #create dict
        d = {}
        for v in Arrangements(set,k):
            tmp_dict = {}
            #add edges of dimension i
            for i in xrange(1,k):
                #swap 0th and ith element
                v[0], v[i] = v[i], v[0]
                #convert to str and add to list
                vert = "".join(v)
                tmp_dict[vert] = None
                #swap back
                v[0], v[i] = v[i], v[0]
            #add other edges
            tmp_bit = v[0]
            for i in set:
                #check if external
                if not (i in v):
                    v[0] = i
                    #add edge
                    vert = "".join(v)
                    tmp_dict[vert] = None
                v[0] = tmp_bit
            d["".join(v)] = tmp_dict
        return graph.Graph(d, name="(%d,%d)-star"%(n,k))

    def NStarGraph(self,n):
        r"""
        Returns the n-star graph.

        The vertices of the n-star graph are the set of permutations on n
        symbols. There is an edge between two vertices if their labels differ
        only in the first and one other position.

        INPUT:

        -  ``n``

        EXAMPLES::

            sage: g = graphs.NStarGraph(4)
            sage: g.plot() # long time

        REFERENCES:

        - S.B. Akers, D. Horel and B. Krishnamurthy, The star graph: An
          attractive alternative to the previous n-cube. In: Proc. Internat.
          Conf. on Parallel Processing (1987), pp. 393--400.

        AUTHORS:

        - Michael Yurko (2009-09-01)
        """
        from sage.combinat.permutation import Permutations
        #set from which to permute
        set = [str(i) for i in xrange(1,n+1)]
        #create dictionary of lists
        #vertices are adjacent if the first element
        #is swapped with the ith element
        d = {}
        for v in Permutations(set):
            tmp_dict = {}
            for i in xrange(1,n):
                if v[0] != v[i]:
                    #swap 0th and ith element
                    v[0], v[i] = v[i], v[0]
                    #convert to str and add to list
                    vert = "".join(v)
                    tmp_dict[vert] = None
                    #swap back
                    v[0], v[i] = v[i], v[0]
            d["".join(v)] = tmp_dict
        return graph.Graph(d, name = "%d-star"%n)

    def OddGraph(self,n):
        r"""
        Returns the Odd Graph with parameter `n`.

        The Odd Graph with parameter `n` is defined as the
        Kneser Graph with parameters `2n-1,n-1`.
        Equivalently, the Odd Graph is the graph whose vertices
        are the `n-1`-subsets of `[0,1,\dots,2(n-1)]`, and such
        that two vertices are adjacent if their corresponding sets
        are disjoint.

        For example, the Petersen Graph can be defined
        as the Odd Graph with parameter `3`.

        EXAMPLE::

            sage: OG=graphs.OddGraph(3)
            sage: print OG.vertices()
            [{4, 5}, {1, 3}, {2, 5}, {2, 3}, {3, 4}, {3, 5}, {1, 4}, {1, 5}, {1, 2}, {2, 4}]
            sage: P=graphs.PetersenGraph()
            sage: P.is_isomorphic(OG)
            True

        TESTS::

            sage: KG=graphs.OddGraph(1)
            Traceback (most recent call last):
            ...
            ValueError: Parameter n should be an integer strictly greater than 1
        """

        if not n>1:
            raise ValueError, "Parameter n should be an integer strictly greater than 1"
        g = self.KneserGraph(2*n-1,n-1)
        g.name("Odd Graph with parameter %s" % n)
        return g

    def PaleyGraph(self,q):
        r"""
        Paley graph with `q` vertices

        Parameter `q` must be the power of a prime number and congruent
        to 1 mod 4.

        EXAMPLES::

            sage: G=graphs.PaleyGraph(9);G
            Paley graph with parameter 9: Graph on 9 vertices
            sage: G.is_regular()
            True

        A Paley graph is always self-complementary::

            sage: G.complement().is_isomorphic(G)
            True
        """
        from sage.rings.finite_rings.integer_mod import mod
        from sage.rings.finite_rings.constructor import FiniteField
        assert q.is_prime_power(), "Parameter q must be a prime power"
        assert mod(q,4)==1, "Parameter q must be congruent to 1 mod 4"
        g = graph.Graph([FiniteField(q,'a'), lambda i,j: (i-j).is_square()],
        loops=False, name = "Paley graph with parameter %d"%q)
        return g


    def PermutationGraph(self, second_permutation, first_permutation = None):
        r"""
        Builds a permutation graph from one (or two) permutations.

        General definition

        A Permutation Graph can be encoded by a permutation `\sigma`
        of `0, ..., n`. It is then built in the following way :

          Take two horizontal lines in the euclidean plane, and mark points `0,
          ..., n` from left to right on the first of them. On the second one,
          still from left to right, mark point in the order in which they appear
          in `\sigma`. Now, link by a segment the two points marked with 1, then
          link together the points marked with 2, and so on. The permutation
          graph defined by the permutation is the intersection graph of those
          segments : there exists a point in this graph for each element from
          `1` to `n`, two vertices `i, j` being adjacent if the segments `i` and
          `j` cross each other.

        The set of edges of the resulting graph is equal to the set of
        inversions of the inverse of the given permutation.

        INPUT:

        - ``second_permutation`` -- the permutation from which the graph should
          be built. It corresponds to the ordering of the elements on the second
          line (see previous definition)

        - ``first_permutation`` (optional) -- the ordering of the elements on
          the *first* line. This is useful when the elements have no natural
          ordering, for instance when they are strings, or tuples, or anything
          else.

          When ``first_permutation == None`` (default), it is set to be equal to
          ``sorted(second_permutation)``, which just yields the expected
          ordering when the elements of the graph are integers.

        .. SEEALSO:

          - Recognition of Permutation graphs in the :mod:`comparability module
            <sage.graphs.comparability>`.

          - Drawings of permutation graphs as intersection graphs of segments is
            possible through the
            :meth:`~sage.combinat.permutation.Permutation_class.show` method of
            :class:`~sage.combinat.permutation.Permutation` objects.

            The correct argument to use in this case is ``show(representation =
            "braid")``.

          - :meth:`~sage.combinat.permutation.Permutation_class.inversions`

        EXAMPLE::

            sage: p = Permutations(5).random_element()
            sage: edges = graphs.PermutationGraph(p).edges(labels =False)
            sage: set(edges) == set(map(lambda (x,y) : (x+1,y+1),p.inverse().inversions()))
            True

        TESTS::

            sage: graphs.PermutationGraph([1, 2, 3], [4, 5, 6])
            Traceback (most recent call last):
            ...
            ValueError: The two permutations do not contain the same set of elements ...
        """
        if first_permutation == None:
            first_permutation = sorted(second_permutation)
        else:
            if set(second_permutation) != set(first_permutation):
                raise ValueError("The two permutations do not contain the same "+
                                 "set of elements ! It is going to be pretty "+
                                 "hard to define a permutation graph from that !")

        vertex_to_index = {}
        for i, v in enumerate(first_permutation):
            vertex_to_index[v] = i+1

        from sage.combinat.permutation import Permutation
        p2 = Permutation(map(lambda x:vertex_to_index[x], second_permutation))
        p1 = Permutation(map(lambda x:vertex_to_index[x], first_permutation))
        p2 = p2 * p1.inverse()
        p2 = p2.inverse()

        g = graph.Graph(name="Permutation graph for "+str(second_permutation))
        g.add_vertices(second_permutation)

        for u,v in p2.inversions():
            g.add_edge(first_permutation[u], first_permutation[v])

        return g

################################################################################
#   Pseudofractal Graphs
################################################################################

    def DorogovtsevGoltsevMendesGraph(self, n):
        """
        Construct the n-th generation of the Dorogovtsev-Goltsev-Mendes
        graph.

        EXAMPLE::

            sage: G = graphs.DorogovtsevGoltsevMendesGraph(8)
            sage: G.size()
            6561

        REFERENCE:

        - [1] Dorogovtsev, S. N., Goltsev, A. V., and Mendes, J.
          F. F., Pseudofractal scale-free web, Phys. Rev. E 066122
          (2002).
        """
        import networkx
        return graph.Graph(networkx.dorogovtsev_goltsev_mendes_graph(n),\
               name="Dorogovtsev-Goltsev-Mendes Graph, %d-th generation"%n)

################################################################################
#   Random Graphs
################################################################################

    def RandomGNP(self, n, p, seed=None, fast=True, method='Sage'):
        r"""
        Returns a random graph on `n` nodes. Each edge is inserted independently
        with probability `p`.

        INPUTS:

        - ``n`` -- number of nodes of the digraph

        - ``p`` -- probability of an edge

        - ``seed`` -- integer seed for random number generator (default=None).

        - ``fast`` -- boolean set to True (default) to use the algorithm with
          time complexity in `O(n+m)` proposed in [3]_. It is designed for
          generating large sparse graphs. It is faster than other methods for
          *LARGE* instances (try it to know whether it is useful for you).

        - ``method`` -- By default (```method='Sage'``), this function uses the
          method implemented in ```sage.graphs.graph_generators_pyx.pyx``. When
          ``method='networkx'``, this function calls the NetworkX function
          ``fast_gnp_random_graph``, unless ``fast=False``, then
          ``gnp_random_graph``. Try them to know which method is the best for
          you. The ``fast`` parameter is not taken into account by the 'Sage'
          method so far.

        REFERENCES:

        .. [1] P. Erdos and A. Renyi. On Random Graphs, Publ.  Math. 6, 290 (1959).

        .. [2] E. N. Gilbert. Random Graphs, Ann. Math.  Stat., 30, 1141 (1959).

        .. [3] V. Batagelj and U. Brandes. Efficient generation of large
               random networks. Phys. Rev. E, 71, 036113, 2005.

        PLOTTING: When plotting, this graph will use the default spring-layout
        algorithm, unless a position dictionary is specified.

        EXAMPLES: We show the edge list of a random graph on 6 nodes with
        probability `p = .4`::

            sage: set_random_seed(0)
            sage: graphs.RandomGNP(6, .4).edges(labels=False)
            [(0, 1), (0, 5), (1, 2), (2, 4), (3, 4), (3, 5), (4, 5)]

        We plot a random graph on 12 nodes with probability
        `p = .71`::

            sage: gnp = graphs.RandomGNP(12,.71)
            sage: gnp.show() # long time

        We view many random graphs using a graphics array::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.RandomGNP(i+3,.43)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show() # long time
            sage: graphs.RandomGNP(4,1)
            Complete graph: Graph on 4 vertices

        TESTS::

            sage: graphs.RandomGNP(50,.2,method=50)
            Traceback (most recent call last):
            ...
            ValueError: 'method' must be equal to 'networkx' or to 'Sage'.
            sage: set_random_seed(0)
            sage: graphs.RandomGNP(50,.2, method="Sage").size()
            243
            sage: graphs.RandomGNP(50,.2, method="networkx").size()
            258
        """
        if n < 0:
            raise ValueError("The number of nodes must be positive or null.")
        if 0.0 > p or 1.0 < p:
            raise ValueError("The probability p must be in [0..1].")

        if seed is None:
            seed = current_randstate().long_seed()
        if p == 1:
            return graphs.CompleteGraph(n)

        if method == 'networkx':
            import networkx
            if fast:
                G = networkx.fast_gnp_random_graph(n, p, seed=seed)
            else:
                G = networkx.gnp_random_graph(n, p, seed=seed)
            return graph.Graph(G)
        elif method in ['Sage', 'sage']:
            # We use the Sage generator
            from sage.graphs.graph_generators_pyx import RandomGNP as sageGNP
            return sageGNP(n, p)
        else:
            raise ValueError("'method' must be equal to 'networkx' or to 'Sage'.")

    def RandomBarabasiAlbert(self, n, m, seed=None):
        u"""
        Return a random graph created using the Barabasi-Albert preferential
        attachment model.

        A graph with m vertices and no edges is initialized, and a graph of n
        vertices is grown by attaching new vertices each with m edges that are
        attached to existing vertices, preferentially with high degree.

        INPUT:

        - ``n`` - number of vertices in the graph

        - ``m`` - number of edges to attach from each new node

        - ``seed`` - for random number generator

        EXAMPLES:

        We show the edge list of a random graph on 6 nodes with m = 2.

        ::

            sage: graphs.RandomBarabasiAlbert(6,2).edges(labels=False)
            [(0, 2), (0, 3), (0, 4), (1, 2), (2, 3), (2, 4), (2, 5), (3, 5)]

        We plot a random graph on 12 nodes with m = 3.

        ::

            sage: ba = graphs.RandomBarabasiAlbert(12,3)
            sage: ba.show()  # long time

        We view many random graphs using a graphics array::

            sage: g = []
            sage: j = []
            sage: for i in range(1,10):
            ...    k = graphs.RandomBarabasiAlbert(i+3, 3)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show()  # long time

        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return graph.Graph(networkx.barabasi_albert_graph(n,m,seed=seed))

    def RandomBipartite(self, n1,n2, p):
        r"""
        Returns a bipartite graph with `n1+n2` vertices
        such that any edge from `[n1]` to `[n2]` exists
        with probability `p`.

        INPUT:

            - ``n1,n2`` : Cardinalities of the two sets
            - ``p``   : Probability for an edge to exist


        EXAMPLE::

            sage: g=graphs.RandomBipartite(5,2,0.5)
            sage: g.vertices()
            [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (1, 0), (1, 1)]

        TESTS::

            sage: g=graphs.RandomBipartite(5,-3,0.5)
            Traceback (most recent call last):
            ...
            ValueError: n1 and n2 should be integers strictly greater than 0
            sage: g=graphs.RandomBipartite(5,3,1.5)
            Traceback (most recent call last):
            ...
            ValueError: Parameter p is a probability, and so should be a real value between 0 and 1

        Trac ticket #12155::

            sage: graphs.RandomBipartite(5,6,.2).complement()
            complement(Random bipartite graph of size 5+6 with edge probability 0.200000000000000): Graph on 11 vertices
        """
        if not (p>=0 and p<=1):
            raise ValueError, "Parameter p is a probability, and so should be a real value between 0 and 1"
        if not (n1>0 and n2>0):
            raise ValueError, "n1 and n2 should be integers strictly greater than 0"

        from numpy.random import uniform
        from sage.graphs.all import Graph

        g=Graph(name="Random bipartite graph of size "+str(n1) +"+"+str(n2)+" with edge probability "+str(p))

        S1=[(0,i) for i in range(n1)]
        S2=[(1,i) for i in range(n2)]
        g.add_vertices(S1)
        g.add_vertices(S2)

        for w in range(n2):
            for v in range(n1):
                if uniform()<=p :
                    g.add_edge((0,v),(1,w))

        pos = {}
        for i in range(n1):
            pos[(0,i)] = (0, i/(n1-1.0))
        for i in range(n2):
            pos[(1,i)] = (1, i/(n2-1.0))

        g.set_pos(pos)

        return g

    def RandomGNM(self, n, m, dense=False, seed=None):
        """
        Returns a graph randomly picked out of all graphs on n vertices
        with m edges.

        INPUT:


        -  ``n`` - number of vertices.

        -  ``m`` - number of edges.

        -  ``dense`` - whether to use NetworkX's
           dense_gnm_random_graph or gnm_random_graph


        EXAMPLES: We show the edge list of a random graph on 5 nodes with
        10 edges.

        ::

            sage: graphs.RandomGNM(5, 10).edges(labels=False)
            [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]

        We plot a random graph on 12 nodes with m = 12.

        ::

            sage: gnm = graphs.RandomGNM(12, 12)
            sage: gnm.show()  # long time

        We view many random graphs using a graphics array::

            sage: g = []
            sage: j = []
            sage: for i in range(9):
            ...    k = graphs.RandomGNM(i+3, i^2-i)
            ...    g.append(k)
            ...
            sage: for i in range(3):
            ...    n = []
            ...    for m in range(3):
            ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
            ...    j.append(n)
            ...
            sage: G = sage.plot.graphics.GraphicsArray(j)
            sage: G.show()  # long time
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        if dense:
            return graph.Graph(networkx.dense_gnm_random_graph(n, m, seed=seed))
        else:
            return graph.Graph(networkx.gnm_random_graph(n, m, seed=seed))

    def RandomNewmanWattsStrogatz(self, n, k, p, seed=None):
        """
        Returns a Newman-Watts-Strogatz small world random graph on n
        vertices.

        From the NetworkX documentation: First create a ring over n nodes.
        Then each node in the ring is connected with its k nearest
        neighbors. Then shortcuts are created by adding new edges as
        follows: for each edge u-v in the underlying "n-ring with k nearest
        neighbors"; with probability p add a new edge u-w with
        randomly-chosen existing node w. In contrast with
        watts_strogatz_graph(), no edges are removed.

        INPUT:


        -  ``n`` - number of vertices.

        -  ``k`` - each vertex is connected to its k nearest
           neighbors

        -  ``p`` - the probability of adding a new edge for
           each edge

        -  ``seed`` - for the random number generator


        EXAMPLE: We show the edge list of a random graph on 7 nodes with 2
        "nearest neighbors" and probability `p = 0.2`::

            sage: graphs.RandomNewmanWattsStrogatz(7, 2, 0.2).edges(labels=False)
            [(0, 1), (0, 2), (0, 3), (0, 6), (1, 2), (2, 3), (2, 4), (3, 4), (3, 6), (4, 5), (5, 6)]

        ::

            sage: G = graphs.RandomNewmanWattsStrogatz(12, 2, .3)
            sage: G.show()  # long time

        REFERENCE:

        - [1] Newman, M.E.J., Watts, D.J. and Strogatz, S.H.  Random
          graph models of social networks. Proc. Nat. Acad. Sci. USA
          99, 2566-2572.
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return graph.Graph(networkx.newman_watts_strogatz_graph(n, k, p, seed=seed))

    def RandomHolmeKim(self, n, m, p, seed=None):
        """
        Returns a random graph generated by the Holme and Kim algorithm for
        graphs with power law degree distribution and approximate average
        clustering.

        INPUT:


        -  ``n`` - number of vertices.

        -  ``m`` - number of random edges to add for each new
           node.

        -  ``p`` - probability of adding a triangle after
           adding a random edge.

        -  ``seed`` - for the random number generator.


        From the NetworkX documentation: The average clustering has a hard
        time getting above a certain cutoff that depends on m. This cutoff
        is often quite low. Note that the transitivity (fraction of
        triangles to possible triangles) seems to go down with network
        size. It is essentially the Barabasi-Albert growth model with an
        extra step that each random edge is followed by a chance of making
        an edge to one of its neighbors too (and thus a triangle). This
        algorithm improves on B-A in the sense that it enables a higher
        average clustering to be attained if desired. It seems possible to
        have a disconnected graph with this algorithm since the initial m
        nodes may not be all linked to a new node on the first iteration
        like the BA model.

        EXAMPLE: We show the edge list of a random graph on 8 nodes with 2
        random edges per node and a probability `p = 0.5` of
        forming triangles.

        ::

            sage: graphs.RandomHolmeKim(8, 2, 0.5).edges(labels=False)
            [(0, 2), (0, 5), (1, 2), (1, 3), (2, 3), (2, 4), (2, 6), (2, 7),
             (3, 4), (3, 6), (3, 7), (4, 5)]

        ::

            sage: G = graphs.RandomHolmeKim(12, 3, .3)
            sage: G.show()  # long time

        REFERENCE:

        - [1] Holme, P. and Kim, B.J. Growing scale-free networks with
          tunable clustering, Phys. Rev. E (2002). vol 65, no 2,
          026107.
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return graph.Graph(networkx.powerlaw_cluster_graph(n, m, p, seed=seed))

    def RandomInterval(self,n):
        """
        Returns a random interval graph.

        An interval graph is built from a list `(a_i,b_i)_{1\leq i \leq n}`
        of intervals : to each interval of the list is associated one
        vertex, two vertices being adjacent if the two corresponding
        intervals intersect.

        A random interval graph of order `n` is generated by picking
        random values for the `(a_i,b_j)`, each of the two coordinates
        being generated from the uniform distribution on the interval
        `[0,1]`.

        This definitions follows [boucheron2001]_.

        .. NOTE::

            The vertices are named 0, 1, 2, and so on. The intervals
            used to create the graph are saved with the graph and can
            be recovered using ``get_vertex()`` or ``get_vertices()``.

        INPUT:

        - ``n`` (integer) -- the number of vertices in the random
          graph.

        EXAMPLE:

        As for any interval graph, the chromatic number is equal to
        the clique number ::

            sage: g = graphs.RandomInterval(8)
            sage: g.clique_number() == g.chromatic_number()
            True

        REFERENCE:

        .. [boucheron2001] Boucheron, S. and FERNANDEZ de la VEGA, W.,
           On the Independence Number of Random Interval Graphs,
           Combinatorics, Probability and Computing v10, issue 05,
           Pages 385--396,
           Cambridge Univ Press, 2001
        """

        from sage.misc.prandom import random

        intervals = [tuple(sorted((random(), random()))) for i in range(n)]

        return self.IntervalGraph(intervals)

    def IntervalGraph(self,intervals):
        r"""
        Returns the graph corresponding to the given intervals.

        An interval graph is built from a list `(a_i,b_i)_{1\leq i \leq n}`
        of intervals : to each interval of the list is associated one
        vertex, two vertices being adjacent if the two corresponding
        (closed) intervals intersect.

        INPUT:

        - ``intervals`` -- the list of pairs `(a_i,b_i)`
          defining the graph.

        .. NOTE::

            * The vertices are named 0, 1, 2, and so on. The
              intervals used to create the graph are saved with the
              graph and can be recovered using ``get_vertex()`` or
              ``get_vertices()``.

            * The intervals `(a_i,b_i)` need not verify `a_i<b_i`.

        EXAMPLE:

        The following line creates the sequence of intervals
        `(i, i+2)` for i in `[0, ..., 8]`::

            sage: intervals = [(i,i+2) for i in range(9)]

        In the corresponding graph... ::

            sage: g = graphs.IntervalGraph(intervals)
            sage: g.get_vertex(3)
            (3, 5)
            sage: neigh = g.neighbors(3)
            sage: for v in neigh: print g.get_vertex(v)
            (1, 3)
            (2, 4)
            (4, 6)
            (5, 7)

        The is_interval() method verifies that this graph is an interval
        graph. ::

            sage: g.is_interval()
            True

        The intervals in the list need not be distinct. ::

            sage: intervals = [ (1,2), (1,2), (1,2), (2,3), (3,4) ]
            sage: g = graphs.IntervalGraph(intervals)
            sage: g.clique_maximum()
            [0, 1, 2, 3]
            sage: g.get_vertices()
            {0: (1, 2), 1: (1, 2), 2: (1, 2), 3: (2, 3), 4: (3, 4)}

        """

        n = len(intervals)
        g = graph.Graph(n)

        edges = []

        for i in range(n-1):
            I = intervals[i]
            for j in range(i+1,n):
                J = intervals[j]
                if max(I) < min(J) or max(J) < min(I): continue
                edges.append((i,j))

        g.add_edges(edges)

        rep = dict( zip(range(n),intervals) )
        g.set_vertices(rep)

        return g


    def RandomLobster(self, n, p, q, seed=None):
        """
        Returns a random lobster.

        A lobster is a tree that reduces to a caterpillar when pruning all
        leaf vertices. A caterpillar is a tree that reduces to a path when
        pruning all leaf vertices (q=0).

        INPUT:


        -  ``n`` - expected number of vertices in the backbone

        -  ``p`` - probability of adding an edge to the
           backbone

        -  ``q`` - probability of adding an edge (claw) to the
           arms

        -  ``seed`` - for the random number generator


        EXAMPLE: We show the edge list of a random graph with 3 backbone
        nodes and probabilities `p = 0.7` and `q = 0.3`::

            sage: graphs.RandomLobster(3, 0.7, 0.3).edges(labels=False)
            [(0, 1), (1, 2)]

        ::

            sage: G = graphs.RandomLobster(9, .6, .3)
            sage: G.show()  # long time
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return graph.Graph(networkx.random_lobster(n, p, q, seed=seed))

    def RandomTree(self, n):
        """
        Returns a random tree on `n` nodes numbered `0` through `n-1`.

        By Cayley's theorem, there are `n^{n-2}` trees with vertex
        set `\{0,1,...,n-1\}`. This constructor chooses one of these uniformly
        at random.

        ALGORITHM:

        The algoritm works by generating an `(n-2)`-long
        random sequence of numbers chosen independently and uniformly
        from `\{0,1,\ldots,n-1\}` and then applies an inverse
        Prufer transformation.

        INPUT:

        -  ``n`` - number of vertices in the tree

        EXAMPLE::

            sage: G = graphs.RandomTree(10)
            sage: G.is_tree()
            True
            sage: G.show() # long

        TESTS:

        Ensuring that we encounter no unexpected surprise ::

            sage: all( graphs.RandomTree(10).is_tree()
            ...        for i in range(100) )
            True

        """
        from sage.misc.prandom import randint
        g = graph.Graph()

        # create random Prufer code
        code = [ randint(0,n-1) for i in xrange(n-2) ]

        # We count the number of symbols of each type.
        # count[k] is the no. of times k appears in code
        #
        # (count[k] is set to -1 when the corresponding vertex is not
        # available anymore)
        count = [ 0 for i in xrange(n) ]
        for k in code:
            count[k] += 1

        g.add_vertices(range(n))

        for s in code:
            for x in range(n):
                if count[x] == 0:
                    break

            count[x] = -1
            g.add_edge(x,s)
            count[s] -= 1

        # Adding as an edge the last two available vertices
        last_edge = [ v for v in range(n) if count[v] != -1 ]
        g.add_edge(last_edge)

        return g

    def RandomTreePowerlaw(self, n, gamma=3, tries=100, seed=None):
        """
        Returns a tree with a power law degree distribution. Returns False
        on failure.

        From the NetworkX documentation: A trial power law degree sequence
        is chosen and then elements are swapped with new elements from a
        power law distribution until the sequence makes a tree (size = order
        - 1).

        INPUT:


        -  ``n`` - number of vertices

        -  ``gamma`` - exponent of power law

        -  ``tries`` - number of attempts to adjust sequence to
           make a tree

        -  ``seed`` - for the random number generator


        EXAMPLE: We show the edge list of a random graph with 10 nodes and
        a power law exponent of 2.

        ::

            sage: graphs.RandomTreePowerlaw(10, 2).edges(labels=False)
            [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (6, 8), (6, 9)]

        ::

            sage: G = graphs.RandomTreePowerlaw(15, 2)
            sage: if G:
            ...    G.show()  # random output, long time
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        try:
            return graph.Graph(networkx.random_powerlaw_tree(n, gamma, seed=seed, tries=tries))
        except:
            return False

    def RandomRegular(self, d, n, seed=None):
        """
        Returns a random d-regular graph on n vertices, or returns False on
        failure.

        Since every edge is incident to two vertices, n\*d must be even.

        INPUT:


        -  ``n`` - number of vertices

        -  ``d`` - degree

        -  ``seed`` - for the random number generator


        EXAMPLE: We show the edge list of a random graph with 8 nodes each
        of degree 3.

        ::

            sage: graphs.RandomRegular(3, 8).edges(labels=False)
            [(0, 1), (0, 4), (0, 7), (1, 5), (1, 7), (2, 3), (2, 5), (2, 6), (3, 4), (3, 6), (4, 5), (6, 7)]

        ::

            sage: G = graphs.RandomRegular(3, 20)
            sage: if G:
            ...    G.show()  # random output, long time

        REFERENCES:

        - [1] Kim, Jeong Han and Vu, Van H. Generating random regular
          graphs. Proc. 35th ACM Symp. on Thy. of Comp. 2003, pp
          213-222. ACM Press, San Diego, CA, USA.
          http://doi.acm.org/10.1145/780542.780576

        - [2] Steger, A. and Wormald, N. Generating random regular
          graphs quickly. Prob. and Comp. 8 (1999), pp 377-396.
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        try:
            N = networkx.random_regular_graph(d, n, seed=seed)
            if N is False: return False
            return graph.Graph(N, sparse=True)
        except:
            return False

    def RandomShell(self, constructor, seed=None):
        """
        Returns a random shell graph for the constructor given.

        INPUT:


        -  ``constructor`` - a list of 3-tuples (n,m,d), each
           representing a shell

        -  ``n`` - the number of vertices in the shell

        -  ``m`` - the number of edges in the shell

        -  ``d`` - the ratio of inter (next) shell edges to
           intra shell edges

        -  ``seed`` - for the random number generator


        EXAMPLE::

            sage: G = graphs.RandomShell([(10,20,0.8),(20,40,0.8)])
            sage: G.edges(labels=False)
            [(0, 3), (0, 7), (0, 8), (1, 2), (1, 5), (1, 8), (1, 9), (3, 6), (3, 11), (4, 6), (4, 7), (4, 8), (4, 21), (5, 8), (5, 9), (6, 9), (6, 10), (7, 8), (7, 9), (8, 18), (10, 11), (10, 13), (10, 19), (10, 22), (10, 26), (11, 18), (11, 26), (11, 28), (12, 13), (12, 14), (12, 28), (12, 29), (13, 16), (13, 21), (13, 29), (14, 18), (16, 20), (17, 18), (17, 26), (17, 28), (18, 19), (18, 22), (18, 27), (18, 28), (19, 23), (19, 25), (19, 28), (20, 22), (24, 26), (24, 27), (25, 27), (25, 29)]
            sage: G.show()  # long time
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return graph.Graph(networkx.random_shell_graph(constructor, seed=seed))

    def WorldMap(self):
        """
        Returns the Graph of all the countries, in which two countries are adjacent
        in the graph if they have a common boundary.

        This graph has been built from the data available
        in The CIA World Factbook [CIA]_ (2009-08-21).

        The returned graph ``G`` has a member ``G.gps_coordinates``
        equal to a dictionary containing the GPS coordinates
        of each country's capital city.

        EXAMPLE::

            sage: g=graphs.WorldMap()
            sage: g.has_edge("France","Italy")
            True
            sage: g.gps_coordinates["Bolivia"]
            [[17, 'S'], [65, 'W']]
            sage: sorted(g.connected_component_containing_vertex('Ireland'))
            ['Ireland', 'United Kingdom']

        REFERENCE:

        .. [CIA] CIA Factbook 09 https://www.cia.gov/library/publications/the-world-factbook/
        """
        from sage.graphs.all import Graph
        edges = [
            ('Afghanistan', 'China', None), ('Afghanistan', 'Iran', None),
            ('Afghanistan', 'Uzbekistan', None), ('Albania', 'Greece', None),
            ('Albania', 'Kosovo', None), ('Albania', 'Macedonia', None),
            ('Albania', 'Montenegro', None), ('Algeria', 'Morocco', None),
            ('Algeria', 'Tunisia', None), ('Andorra', 'Spain', None),
            ('Angola', 'Democratic Republic of the Congo', None), ('Angola', 'Namibia', None),
            ('Angola', 'Zambia', None), ('Argentina', 'Bolivia', None),
            ('Argentina', 'Brazil', None), ('Argentina', 'Chile', None),
            ('Argentina', 'Paraguay', None), ('Argentina', 'Uruguay', None),
            ('Armenia', 'Georgia', None), ('Armenia', 'Iran', None),
            ('Austria', 'Germany', None), ('Azerbaijan', 'Armenia', None),
            ('Azerbaijan', 'Georgia', None), ('Azerbaijan', 'Iran', None),
            ('Azerbaijan', 'Russia', None), ('Azerbaijan', 'Turkey', None),
            ('Bangladesh', 'Burma', None), ('Belgium', 'Germany', None),
            ('Belgium', 'Netherlands', None), ('Belize', 'Mexico', None),
            ('Benin', 'Burkina Faso', None), ('Benin', 'Niger', None),
            ('Benin', 'Nigeria', None), ('Benin', 'Togo', None),
            ('Bolivia', 'Brazil', None), ('Bolivia', 'Chile', None),
            ('Bolivia', 'Paraguay', None), ('Bolivia', 'Peru', None),
            ('Bosnia and Herzegovina', 'Croatia', None), ('Bosnia and Herzegovina', 'Montenegro', None),
            ('Bosnia and Herzegovina', 'Serbia', None), ('Brazil', 'Colombia', None),
            ('Brazil', 'Guyana', None), ('Brazil', 'Suriname', None),
            ('Brazil', 'Venezuela', None), ('Bulgaria', 'Greece', None),
            ('Bulgaria', 'Macedonia', None), ('Bulgaria', 'Romania', None),
            ('Bulgaria', 'Serbia', None), ('Burkina Faso', 'Mali', None),
            ('Burkina Faso', 'Niger', None), ('Burkina Faso', 'Togo', None),
            ('Burundi', 'Democratic Republic of the Congo', None), ('Cambodia', 'Laos', None),
            ('Cambodia', 'Thailand', None), ('Cambodia', 'Vietnam', None),
            ('Cameroon', 'Central African Republic', None), ('Cameroon', 'Chad', None),
            ('Cameroon', 'Equatorial Guinea', None), ('Cameroon', 'Nigeria', None),
            ('Cameroon', 'Republic of the Congo', None), ('Canada', 'United States', None),
            ('Central African Republic', 'Chad', None), ('Central African Republic', 'Democratic Republic of the Congo', None),
            ('Central African Republic', 'Sudan', None), ('Chad', 'Niger', None),
            ('Chad', 'Nigeria', None), ('Chad', 'Sudan', None),
            ('China', 'Bhutan', None), ('China', 'Burma', None),
            ('China', 'Hong Kong', None), ('China', 'Kazakhstan', None),
            ('China', 'Kyrgyzstan', None), ('China', 'Mongolia', None),
            ('China', 'Nepal', None), ('China', 'North Korea', None),
            ('China', 'Russia', None), ('China', 'Vietnam', None),
            ('Colombia', 'Venezuela', None), ('Costa Rica', 'Nicaragua', None),
            ("Cote d'Ivoire", 'Burkina Faso', None), ("Cote d'Ivoire", 'Guinea', None),
            ("Cote d'Ivoire", 'Mali', None), ('Cyprus', 'Akrotiri', None),
            ('Cyprus', 'Dhekelia', None), ('Czech Republic', 'Austria', None),
            ('Czech Republic', 'Germany', None), ('Czech Republic', 'Poland', None),
            ('Democratic Republic of the Congo', 'Zambia', None), ('Denmark', 'Germany', None),
            ('Djibouti', 'Eritrea', None), ('Dominican Republic', 'Haiti', None),
            ('Ecuador', 'Colombia', None), ('El Salvador', 'Honduras', None),
            ('Ethiopia', 'Djibouti', None), ('Ethiopia', 'Eritrea', None),
            ('Ethiopia', 'Kenya', None), ('Ethiopia', 'Somalia', None),
            ('Ethiopia', 'Sudan', None), ('Finland', 'Russia', None),
            ('Finland', 'Sweden', None), ('France', 'Andorra', None),
            ('France', 'Belgium', None), ('France', 'Brazil', None),
            ('France', 'Germany', None), ('France', 'Italy', None),
            ('France', 'Luxembourg', None), ('France', 'Spain', None),
            ('France', 'Suriname', None), ('France', 'Switzerland', None),
            ('Gabon', 'Cameroon', None), ('Gabon', 'Equatorial Guinea', None),
            ('Gabon', 'Republic of the Congo', None), ('Gaza Strip', 'Egypt', None),
            ('Gaza Strip', 'Israel', None), ('Ghana', 'Burkina Faso', None),
            ('Ghana', "Cote d'Ivoire", None), ('Ghana', 'Togo', None),
            ('Gibraltar', 'Spain', None), ('Guatemala', 'Belize', None),
            ('Guatemala', 'El Salvador', None), ('Guatemala', 'Honduras', None),
            ('Guatemala', 'Mexico', None), ('Guinea', 'Sierra Leone', None),
            ('Guinea-Bissau', 'Guinea', None), ('Guinea-Bissau', 'Senegal', None),
            ('Honduras', 'Nicaragua', None), ('Hungary', 'Austria', None),
            ('Hungary', 'Croatia', None), ('Hungary', 'Serbia', None),
            ('India', 'Bangladesh', None), ('India', 'Bhutan', None),
            ('India', 'Burma', None), ('India', 'China', None),
            ('India', 'Nepal', None), ('Indonesia', 'Papua New Guinea', None),
            ('Iran', 'Iraq', None), ('Ireland', 'United Kingdom', None),
            ('Israel', 'Egypt', None), ('Italy', 'Austria', None),
            ('Jordan', 'Iraq', None), ('Jordan', 'Israel', None),
            ('Jordan', 'Syria', None), ('Jordan', 'West Bank', None),
            ('Kazakhstan', 'Kyrgyzstan', None), ('Kenya', 'Somalia', None),
            ('Kenya', 'Sudan', None), ('Kenya', 'Uganda', None),
            ('Kosovo', 'Macedonia', None), ('Kosovo', 'Serbia', None),
            ('Kuwait', 'Iraq', None), ('Laos', 'Burma', None),
            ('Laos', 'China', None), ('Laos', 'Thailand', None),
            ('Laos', 'Vietnam', None), ('Latvia', 'Belarus', None),
            ('Latvia', 'Estonia', None), ('Lebanon', 'Israel', None),
            ('Lesotho', 'South Africa', None), ('Liberia', "Cote d'Ivoire", None),
            ('Liberia', 'Guinea', None), ('Liberia', 'Sierra Leone', None),
            ('Libya', 'Algeria', None), ('Libya', 'Chad', None),
            ('Libya', 'Egypt', None), ('Libya', 'Niger', None),
            ('Libya', 'Sudan', None), ('Libya', 'Tunisia', None),
            ('Liechtenstein', 'Austria', None), ('Liechtenstein', 'Switzerland', None),
            ('Lithuania', 'Belarus', None), ('Lithuania', 'Latvia', None),
            ('Lithuania', 'Poland', None), ('Lithuania', 'Russia', None),
            ('Luxembourg', 'Belgium', None), ('Luxembourg', 'Germany', None),
            ('Macau', 'China', None), ('Macedonia', 'Greece', None),
            ('Macedonia', 'Serbia', None), ('Malaysia', 'Brunei', None),
            ('Malaysia', 'Indonesia', None), ('Malaysia', 'Thailand', None),
            ('Mali', 'Algeria', None), ('Mali', 'Guinea', None),
            ('Mali', 'Niger', None), ('Mali', 'Senegal', None),
            ('Mauritania', 'Algeria', None), ('Mauritania', 'Mali', None),
            ('Mauritania', 'Senegal', None), ('Mauritania', 'Western Sahara', None),
            ('Monaco', 'France', None), ('Montenegro', 'Croatia', None),
            ('Montenegro', 'Kosovo', None), ('Montenegro', 'Serbia', None),
            ('Morocco', 'Spain', None), ('Mozambique', 'Malawi', None),
            ('Mozambique', 'Zambia', None), ('Mozambique', 'Zimbabwe', None),
            ('Namibia', 'Botswana', None), ('Namibia', 'Zambia', None),
            ('Netherlands', 'Germany', None), ('Niger', 'Algeria', None),
            ('Niger', 'Nigeria', None), ('Norway', 'Finland', None),
            ('Norway', 'Russia', None), ('Norway', 'Sweden', None),
            ('Oman', 'United Arab Emirates', None), ('Oman', 'Yemen', None),
            ('Pakistan', 'Afghanistan', None), ('Pakistan', 'China', None),
            ('Pakistan', 'India', None), ('Pakistan', 'Iran', None),
            ('Panama', 'Colombia', None), ('Panama', 'Costa Rica', None),
            ('Paraguay', 'Brazil', None), ('Peru', 'Brazil', None),
            ('Peru', 'Chile', None), ('Peru', 'Colombia', None),
            ('Peru', 'Ecuador', None), ('Poland', 'Belarus', None),
            ('Poland', 'Germany', None), ('Portugal', 'Spain', None),
            ('Republic of the Congo', 'Angola', None), ('Republic of the Congo', 'Central African Republic', None),
            ('Republic of the Congo', 'Democratic Republic of the Congo', None), ('Romania', 'Hungary', None),
            ('Romania', 'Moldova', None), ('Romania', 'Serbia', None),
            ('Russia', 'Belarus', None), ('Russia', 'Estonia', None),
            ('Russia', 'Georgia', None), ('Russia', 'Kazakhstan', None),
            ('Russia', 'Latvia', None), ('Russia', 'Mongolia', None),
            ('Russia', 'North Korea', None), ('Russia', 'Poland', None),
            ('Rwanda', 'Burundi', None), ('Rwanda', 'Democratic Republic of the Congo', None),
            ('Rwanda', 'Uganda', None), ('Saint Martin', 'Netherlands Antilles', None),
            ('San Marino', 'Italy', None), ('Saudi Arabia', 'Iraq', None),
            ('Saudi Arabia', 'Jordan', None), ('Saudi Arabia', 'Kuwait', None),
            ('Saudi Arabia', 'Oman', None), ('Saudi Arabia', 'Qatar', None),
            ('Saudi Arabia', 'United Arab Emirates', None), ('Saudi Arabia', 'Yemen', None),
            ('Senegal', 'Guinea', None), ('Serbia', 'Croatia', None),
            ('Slovakia', 'Austria', None), ('Slovakia', 'Czech Republic', None),
            ('Slovakia', 'Hungary', None), ('Slovakia', 'Poland', None),
            ('Slovakia', 'Ukraine', None), ('Slovenia', 'Austria', None),
            ('Slovenia', 'Croatia', None), ('Slovenia', 'Hungary', None),
            ('Slovenia', 'Italy', None), ('Somalia', 'Djibouti', None),
            ('South Africa', 'Botswana', None), ('South Africa', 'Mozambique', None),
            ('South Africa', 'Namibia', None), ('South Africa', 'Zimbabwe', None),
            ('South Korea', 'North Korea', None), ('Sudan', 'Democratic Republic of the Congo', None),
            ('Sudan', 'Egypt', None), ('Sudan', 'Eritrea', None),
            ('Suriname', 'Guyana', None), ('Swaziland', 'Mozambique', None),
            ('Swaziland', 'South Africa', None), ('Switzerland', 'Austria', None),
            ('Switzerland', 'Germany', None), ('Switzerland', 'Italy', None),
            ('Syria', 'Iraq', None), ('Syria', 'Israel', None),
            ('Syria', 'Lebanon', None), ('Tajikistan', 'Afghanistan', None),
            ('Tajikistan', 'China', None), ('Tajikistan', 'Kyrgyzstan', None),
            ('Tajikistan', 'Uzbekistan', None), ('Tanzania', 'Burundi', None),
            ('Tanzania', 'Democratic Republic of the Congo', None), ('Tanzania', 'Kenya', None),
            ('Tanzania', 'Malawi', None), ('Tanzania', 'Mozambique', None),
            ('Tanzania', 'Rwanda', None), ('Tanzania', 'Uganda', None),
            ('Tanzania', 'Zambia', None), ('Thailand', 'Burma', None),
            ('The Gambia', 'Senegal', None), ('Timor-Leste', 'Indonesia', None),
            ('Turkey', 'Armenia', None), ('Turkey', 'Bulgaria', None),
            ('Turkey', 'Georgia', None), ('Turkey', 'Greece', None),
            ('Turkey', 'Iran', None), ('Turkey', 'Iraq', None),
            ('Turkey', 'Syria', None), ('Turkmenistan', 'Afghanistan', None),
            ('Turkmenistan', 'Iran', None), ('Turkmenistan', 'Kazakhstan', None),
            ('Turkmenistan', 'Uzbekistan', None), ('Uganda', 'Democratic Republic of the Congo', None),
            ('Uganda', 'Sudan', None), ('Ukraine', 'Belarus', None),
            ('Ukraine', 'Hungary', None), ('Ukraine', 'Moldova', None),
            ('Ukraine', 'Poland', None), ('Ukraine', 'Romania', None),
            ('Ukraine', 'Russia', None), ('United States', 'Mexico', None),
            ('Uruguay', 'Brazil', None), ('Uzbekistan', 'Kazakhstan', None),
            ('Uzbekistan', 'Kyrgyzstan', None), ('Vatican City', 'Italy', None),
            ('Venezuela', 'Guyana', None), ('West Bank', 'Israel', None),
            ('Western Sahara', 'Algeria', None), ('Western Sahara', 'Morocco', None),
            ('Zambia', 'Malawi', None), ('Zambia', 'Zimbabwe', None),
            ('Zimbabwe', 'Botswana', None)
            ]
        gps_coordinates = {
            'Canada': [[60, 'N'], [95, 'W']],
            'Saint Martin': [[18, 'N'], [63, 'W']],
            'Sao Tome and Principe': [[1, 'N'], [7, 'E']],
            'Turkmenistan': [[40, 'N'], [60, 'E']],
            'Saint Helena': [[15, 'S'], [5, 'W']],
            'Lithuania': [[56, 'N'], [24, 'E']],
            'Cambodia': [[13, 'N'], [105, 'E']],
            'Saint Kitts and Nevis': [[17, 'N'], [62, 'W']],
            'Ethiopia': [[8, 'N'], [38, 'E']],
            'The Gambia': [[13, 'N'], [16, 'W']],
            'Aruba': [[12, 'N'], [69, 'W']],
            'Swaziland': [[26, 'S'], [31, 'E']],
            'Guinea-Bissau': [[12, 'N'], [15, 'W']],
            'Argentina': [[34, 'S'], [64, 'W']],
            'Bolivia': [[17, 'S'], [65, 'W']],
            'Bahamas, The': [[24, 'N'], [76, 'W']],
            'Spratly Islands': [[8, 'N'], [111, 'E']],
            'Ghana': [[8, 'N'], [2, 'W']],
            'Saudi Arabia': [[25, 'N'], [45, 'E']],
            'American Samoa': [[14, 'S'], [170, 'W']],
            'Cocos (Keeling) Islands': [[12, 'S'], [96, 'E']],
            'Slovenia': [[46, 'N'], [14, 'E']],
            'Guatemala': [[15, 'N'], [90, 'W']],
            'Bosnia and Herzegovina': [[44, 'N'], [18, 'E']],
            'Kuwait': [[29, 'N'], [45, 'E']],
            'Jordan': [[31, 'N'], [36, 'E']],
            'Saint Barthelemy': [[17, 'N'], [62, 'W']],
            'Ashmore and Cartier Islands': [[12, 'S'], [123, 'E']],
            'Dominica': [[15, 'N'], [61, 'W']],
            'Liberia': [[6, 'N'], [9, 'W']],
            'Maldives': [[3, 'N'], [73, 'E']],
            'Micronesia, Federated States of': [[6, 'N'], [158, 'E']],
            'Pakistan': [[30, 'N'], [70, 'E']],
            'Oman': [[21, 'N'], [57, 'E']],
            'Tanzania': [[6, 'S'], [35, 'E']],
            'Albania': [[41, 'N'], [20, 'E']],
            'Gabon': [[1, 'S'], [11, 'E']],
            'Niue': [[19, 'S'], [169, 'W']],
            'Monaco': [[43, 'N'], [7, 'E']],
            'Wallis and Futuna': [[13, 'S'], [176, 'W']],
            'New Zealand': [[41, 'S'], [174, 'E']],
            'Yemen': [[15, 'N'], [48, 'E']],
            'Jersey': [[49, 'N'], [2, 'W']],
            'Jamaica': [[18, 'N'], [77, 'W']],
            'Greenland': [[72, 'N'], [40, 'W']],
            'West Bank': [[32, 'N'], [35, 'E']],
            'Macau': [[22, 'N'], [113, 'E']],
            'Jan Mayen': [[71, 'N'], [8, 'W']],
            'United Arab Emirates': [[24, 'N'], [54, 'E']],
            'Guam': [[13, 'N'], [144, 'E']],
            'Uruguay': [[33, 'S'], [56, 'W']],
            'India': [[20, 'N'], [77, 'E']],
            'Azerbaijan': [[40, 'N'], [47, 'E']],
            'Lesotho': [[29, 'S'], [28, 'E']],
            'Saint Vincent and the Grenadines': [[13, 'N'], [61, 'W']],
            'Kenya': [[1, 'N'], [38, 'E']],
            'South Korea': [[37, 'N'], [127, 'E']],
            'Tajikistan': [[39, 'N'], [71, 'E']],
            'Turkey': [[39, 'N'], [35, 'E']],
            'Afghanistan': [[33, 'N'], [65, 'E']],
            'Paraguay': [[23, 'S'], [58, 'W']],
            'Bangladesh': [[24, 'N'], [90, 'E']],
            'Mauritania': [[20, 'N'], [12, 'W']],
            'Solomon Islands': [[8, 'S'], [159, 'E']],
            'Saint Pierre and Miquelon': [[46, 'N'], [56, 'W']],
            'Gaza Strip': [[31, 'N'], [34, 'E']],
            'San Marino': [[43, 'N'], [12, 'E']],
            'French Polynesia': [[15, 'S'], [140, 'W']],
            'France': [[46, 'N'], [2, 'E']],
            'Fiji': [[18, 'S'], [175, 'E']],
            'Rwanda': [[2, 'S'], [30, 'E']],
            'Slovakia': [[48, 'N'], [19, 'E']],
            'Somalia': [[10, 'N'], [49, 'E']],
            'Peru': [[10, 'S'], [76, 'W']],
            'Laos': [[18, 'N'], [105, 'E']],
            'Nauru': [[0, 'S'], [166, 'E']],
            'Seychelles': [[4, 'S'], [55, 'E']],
            'Norway': [[62, 'N'], [10, 'E']],
            "Cote d'Ivoire": [[8, 'N'], [5, 'W']],
            'Cook Islands': [[21, 'S'], [159, 'W']],
            'Benin': [[9, 'N'], [2, 'E']],
            'Western Sahara': [[24, 'N'], [13, 'W']],
            'Cuba': [[21, 'N'], [80, 'W']],
            'Cameroon': [[6, 'N'], [12, 'E']],
            'Montenegro': [[42, 'N'], [19, 'E']],
            'Republic of the Congo': [[1, 'S'], [15, 'E']],
            'Burkina Faso': [[13, 'N'], [2, 'W']],
            'Togo': [[8, 'N'], [1, 'E']],
            'Virgin Islands': [[18, 'N'], [64, 'W']],
            'China': [[35, 'N'], [105, 'E']],
            'Armenia': [[40, 'N'], [45, 'E']],
            'Timor-Leste': [[8, 'S'], [125, 'E']],
            'Dominican Republic': [[19, 'N'], [70, 'W']],
            'Ukraine': [[49, 'N'], [32, 'E']],
            'Bahrain': [[26, 'N'], [50, 'E']],
            'Tonga': [[20, 'S'], [175, 'W']],
            'Finland': [[64, 'N'], [26, 'E']],
            'Libya': [[25, 'N'], [17, 'E']],
            'Cayman Islands': [[19, 'N'], [80, 'W']],
            'Central African Republic': [[7, 'N'], [21, 'E']],
            'New Caledonia': [[21, 'S'], [165, 'E']],
            'Mauritius': [[20, 'S'], [57, 'E']],
            'Liechtenstein': [[47, 'N'], [9, 'E']],
            'Vietnam': [[16, 'N'], [107, 'E']],
            'British Virgin Islands': [[18, 'N'], [64, 'W']],
            'Mali': [[17, 'N'], [4, 'W']],
            'Vatican City': [[41, 'N'], [12, 'E']],
            'Russia': [[60, 'N'], [100, 'E']],
            'Bulgaria': [[43, 'N'], [25, 'E']],
            'United States': [[38, 'N'], [97, 'W']],
            'Romania': [[46, 'N'], [25, 'E']],
            'Angola': [[12, 'S'], [18, 'E']],
            'Chad': [[15, 'N'], [19, 'E']],
            'South Africa': [[29, 'S'], [24, 'E']],
            'Tokelau': [[9, 'S'], [172, 'W']],
            'Turks and Caicos Islands': [[21, 'N'], [71, 'W']],
            'South Georgia and the South Sandwich Islands': [[54, 'S'], [37, 'W']],
            'Sweden': [[62, 'N'], [15, 'E']],
            'Qatar': [[25, 'N'], [51, 'E']],
            'Malaysia': [[2, 'N'], [112, 'E']],
            'Senegal': [[14, 'N'], [14, 'W']],
            'Latvia': [[57, 'N'], [25, 'E']],
            'Clipperton Island': [[10, 'N'], [109, 'W']],
            'Uganda': [[1, 'N'], [32, 'E']],
            'Japan': [[36, 'N'], [138, 'E']],
            'Niger': [[16, 'N'], [8, 'E']],
            'Brazil': [[10, 'S'], [55, 'W']],
            'Faroe Islands': [[62, 'N'], [7, 'W']],
            'Guinea': [[11, 'N'], [10, 'W']],
            'Panama': [[9, 'N'], [80, 'W']],
            'Costa Rica': [[10, 'N'], [84, 'W']],
            'Luxembourg': [[49, 'N'], [6, 'E']],
            'Cape Verde': [[16, 'N'], [24, 'W']],
            'Andorra': [[42, 'N'], [1, 'E']],
            'Gibraltar': [[36, 'N'], [5, 'W']],
            'Ireland': [[53, 'N'], [8, 'W']],
            'Syria': [[35, 'N'], [38, 'E']],
            'Palau': [[7, 'N'], [134, 'E']],
            'Nigeria': [[10, 'N'], [8, 'E']],
            'Ecuador': [[2, 'S'], [77, 'W']],
            'Northern Mariana Islands': [[15, 'N'], [145, 'E']],
            'Brunei': [[4, 'N'], [114, 'E']],
            'Mozambique': [[18, 'S'], [35, 'E']],
            'Australia': [[27, 'S'], [133, 'E']],
            'Iran': [[32, 'N'], [53, 'E']],
            'Algeria': [[28, 'N'], [3, 'E']],
            'Svalbard': [[78, 'N'], [20, 'E']],
            'El Salvador': [[13, 'N'], [88, 'W']],
            'Tuvalu': [[8, 'S'], [178, 'E']],
            'Pitcairn Islands': [[25, 'S'], [130, 'W']],
            'Czech Republic': [[49, 'N'], [15, 'E']],
            'Marshall Islands': [[9, 'N'], [168, 'E']],
            'Chile': [[30, 'S'], [71, 'W']],
            'Puerto Rico': [[18, 'N'], [66, 'W']],
            'Belgium': [[50, 'N'], [4, 'E']],
            'Kiribati': [[1, 'N'], [173, 'E']],
            'Haiti': [[19, 'N'], [72, 'W']],
            'Belize': [[17, 'N'], [88, 'W']],
            'Hong Kong': [[22, 'N'], [114, 'E']],
            'Saint Lucia': [[13, 'N'], [60, 'W']],
            'Georgia': [[42, 'N'], [43, 'E']],
            'Mexico': [[23, 'N'], [102, 'W']],
            'Denmark': [[56, 'N'], [10, 'E']],
            'Poland': [[52, 'N'], [20, 'E']],
            'Moldova': [[47, 'N'], [29, 'E']],
            'Morocco': [[32, 'N'], [5, 'W']],
            'Namibia': [[22, 'S'], [17, 'E']],
            'Mongolia': [[46, 'N'], [105, 'E']],
            'Guernsey': [[49, 'N'], [2, 'W']],
            'Thailand': [[15, 'N'], [100, 'E']],
            'Switzerland': [[47, 'N'], [8, 'E']],
            'Grenada': [[12, 'N'], [61, 'W']],
            'Navassa Island': [[18, 'N'], [75, 'W']],
            'Isle of Man': [[54, 'N'], [4, 'W']],
            'Portugal': [[39, 'N'], [8, 'W']],
            'Estonia': [[59, 'N'], [26, 'E']],
            'Kosovo': [[42, 'N'], [21, 'E']],
            'Norfolk Island': [[29, 'S'], [167, 'E']],
            'Bouvet Island': [[54, 'S'], [3, 'E']],
            'Lebanon': [[33, 'N'], [35, 'E']],
            'Sierra Leone': [[8, 'N'], [11, 'W']],
            'Uzbekistan': [[41, 'N'], [64, 'E']],
            'Tunisia': [[34, 'N'], [9, 'E']],
            'Djibouti': [[11, 'N'], [43, 'E']],
            'Heard Island and McDonald Islands': [[53, 'S'], [72, 'E']],
            'Antigua and Barbuda': [[17, 'N'], [61, 'W']],
            'Spain': [[40, 'N'], [4, 'W']],
            'Colombia': [[4, 'N'], [72, 'W']],
            'Burundi': [[3, 'S'], [30, 'E']],
            'Taiwan': [[23, 'N'], [121, 'E']],
            'Cyprus': [[35, 'N'], [33, 'E']],
            'Barbados': [[13, 'N'], [59, 'W']],
            'Falkland Islands (Islas Malvinas)': [[51, 'S'], [59, 'W']],
            'Madagascar': [[20, 'S'], [47, 'E']],
            'Italy': [[42, 'N'], [12, 'E']],
            'Bhutan': [[27, 'N'], [90, 'E']],
            'Sudan': [[15, 'N'], [30, 'E']],
            'Vanuatu': [[16, 'S'], [167, 'E']],
            'Malta': [[35, 'N'], [14, 'E']],
            'Hungary': [[47, 'N'], [20, 'E']],
            'Democratic Republic of the Congo': [[0, 'N'], [25, 'E']],
            'Netherlands': [[52, 'N'], [5, 'E']],
            'Bermuda': [[32, 'N'], [64, 'W']],
            'Suriname': [[4, 'N'], [56, 'W']],
            'Anguilla': [[18, 'N'], [63, 'W']],
            'Venezuela': [[8, 'N'], [66, 'W']],
            'Netherlands Antilles': [[12, 'N'], [69, 'W']],
            'Israel': [[31, 'N'], [34, 'E']],
            'Paracel Islands': [[16, 'N'], [112, 'E']],
            'Wake Island': [[19, 'N'], [166, 'E']],
            'Indonesia': [[5, 'S'], [120, 'E']],
            'Iceland': [[65, 'N'], [18, 'W']],
            'Zambia': [[15, 'S'], [30, 'E']],
            'Samoa': [[13, 'S'], [172, 'W']],
            'Austria': [[47, 'N'], [13, 'E']],
            'Papua New Guinea': [[6, 'S'], [147, 'E']],
            'Malawi': [[13, 'S'], [34, 'E']],
            'Zimbabwe': [[20, 'S'], [30, 'E']],
            'Germany': [[51, 'N'], [9, 'E']],
            'Dhekelia': [[34, 'N'], [33, 'E']],
            'Kazakhstan': [[48, 'N'], [68, 'E']],
            'Philippines': [[13, 'N'], [122, 'E']],
            'Eritrea': [[15, 'N'], [39, 'E']],
            'Kyrgyzstan': [[41, 'N'], [75, 'E']],
            'Mayotte': [[12, 'S'], [45, 'E']],
            'Iraq': [[33, 'N'], [44, 'E']],
            'Montserrat': [[16, 'N'], [62, 'W']],
            'Coral Sea Islands': [[18, 'S'], [152, 'E']],
            'Macedonia': [[41, 'N'], [22, 'E']],
            'British Indian Ocean Territory': [[6, 'S'], [71, 'E']],
            'North Korea': [[40, 'N'], [127, 'E']],
            'Trinidad and Tobago': [[11, 'N'], [61, 'W']],
            'Akrotiri': [[34, 'N'], [32, 'E']],
            'Guyana': [[5, 'N'], [59, 'W']],
            'Belarus': [[53, 'N'], [28, 'E']],
            'Nepal': [[28, 'N'], [84, 'E']],
            'Burma': [[22, 'N'], [98, 'E']],
            'Honduras': [[15, 'N'], [86, 'W']],
            'Equatorial Guinea': [[2, 'N'], [10, 'E']],
            'Egypt': [[27, 'N'], [30, 'E']],
            'Nicaragua': [[13, 'N'], [85, 'W']],
            'Singapore': [[1, 'N'], [103, 'E']],
            'Serbia': [[44, 'N'], [21, 'E']],
            'Botswana': [[22, 'S'], [24, 'E']],
            'United Kingdom': [[54, 'N'], [2, 'W']],
            'Antarctica': [[90, 'S'], [0, 'E']],
            'Christmas Island': [[10, 'S'], [105, 'E']],
            'Greece': [[39, 'N'], [22, 'E']],
            'Sri Lanka': [[7, 'N'], [81, 'E']],
            'Croatia': [[45, 'N'], [15, 'E']],
            'Comoros': [[12, 'S'], [44, 'E']]
            }
        g = Graph()
        g.add_edges(edges)
        g.gps_coordinates = gps_coordinates
        g.name("World Map")
        return g



################################################################################
#   Graphs with a given degree sequence
################################################################################

    def DegreeSequence(self, deg_sequence):
        """
        Returns a graph with the given degree sequence. Raises a NetworkX
        error if the proposed degree sequence cannot be that of a graph.

        Graph returned is the one returned by the Havel-Hakimi algorithm,
        which constructs a simple graph by connecting vertices of highest
        degree to other vertices of highest degree, resorting the remaining
        vertices by degree and repeating the process. See Theorem 1.4 in
        [1].

        INPUT:


        -  ``deg_sequence`` - a list of integers with each
           entry corresponding to the degree of a different vertex.


        EXAMPLES::

            sage: G = graphs.DegreeSequence([3,3,3,3])
            sage: G.edges(labels=False)
            [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            sage: G.show()  # long time

        ::

            sage: G = graphs.DegreeSequence([3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3])
            sage: G.show()  # long time

        ::

            sage: G = graphs.DegreeSequence([4,4,4,4,4,4,4,4])
            sage: G.show()  # long time

        ::

            sage: G = graphs.DegreeSequence([1,2,3,4,3,4,3,2,3,2,1])
            sage: G.show()  # long time

        REFERENCE:

        - [1] Chartrand, G. and Lesniak, L. Graphs and Digraphs.
          Chapman and Hall/CRC, 1996.
        """
        import networkx
        return graph.Graph(networkx.havel_hakimi_graph([int(i) for i in deg_sequence]))

    def DegreeSequenceBipartite(self, s1 ,s2 ):
        r"""
        Returns a bipartite graph whose two sets have the given
        degree sequences.

        Given two different sequences of degrees `s_1` and `s_2`,
        this functions returns ( if possible ) a bipartite graph
        on sets `A` and `B` such that the vertices in `A` have
        `s_1` as their degree sequence, while `s_2` is the degree
        sequence of the vertices in `B`.

        INPUT:

        - ``s_1`` -- list of integers corresponding to the degree
          sequence of the first set.
        - ``s_2`` -- list of integers corresponding to the degree
          sequence of the second set.

        ALGORITHM:

        This function works through the computation of the matrix
        given by the Gale-Ryser theorem, which is in this case
        the adjacency matrix of the bipartite graph.

        EXAMPLES:

        If we are given as sequences ``[2,2,2,2,2]`` and ``[5,5]``
        we are given as expected the complete bipartite
        graph `K_{2,5}` ::

            sage: g = graphs.DegreeSequenceBipartite([2,2,2,2,2],[5,5])
            sage: g.is_isomorphic(graphs.CompleteBipartiteGraph(5,2))
            True

        Some sequences being incompatible if, for example, their sums
        are different, the functions raises a ``ValueError`` when no
        graph corresponding to the degree sequences exists. ::

            sage: g = graphs.DegreeSequenceBipartite([2,2,2,2,1],[5,5])
            Traceback (most recent call last):
            ...
            ValueError: There exists no bipartite graph corresponding to the given degree sequences

        TESTS:

        Trac ticket #12155::

            sage: graphs.DegreeSequenceBipartite([2,2,2,2,2],[5,5]).complement()
            complement(): Graph on 7 vertices
        """

        from sage.combinat.integer_vector import gale_ryser_theorem
        from sage.graphs.graph import Graph
        from sage.graphs.bipartite_graph import BipartiteGraph

        s1 = sorted(s1, reverse = True)
        s2 = sorted(s2, reverse = True)

        m = gale_ryser_theorem(s1,s2)

        if m is False:
            raise ValueError("There exists no bipartite graph corresponding to the given degree sequences")
        else:
            return Graph(BipartiteGraph(m))

    def DegreeSequenceConfigurationModel(self, deg_sequence, seed=None):
        """
        Returns a random pseudograph with the given degree sequence. Raises
        a NetworkX error if the proposed degree sequence cannot be that of
        a graph with multiple edges and loops.

        One requirement is that the sum of the degrees must be even, since
        every edge must be incident with two vertices.

        INPUT:


        -  ``deg_sequence`` - a list of integers with each
           entry corresponding to the expected degree of a different vertex.

        -  ``seed`` - for the random number generator.


        EXAMPLES::

            sage: G = graphs.DegreeSequenceConfigurationModel([1,1])
            sage: G.adjacency_matrix()
            [0 1]
            [1 0]

        Note: as of this writing, plotting of loops and multiple edges is
        not supported, and the output is allowed to contain both types of
        edges.

        ::

            sage: G = graphs.DegreeSequenceConfigurationModel([3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3])
            sage: G.edges(labels=False)
            [(0, 2), (0, 10), (0, 15), (1, 6), (1, 16), (1, 17), (2, 5), (2, 19), (3, 7), (3, 14), (3, 14), (4, 9), (4, 13), (4, 19), (5, 6), (5, 15), (6, 11), (7, 11), (7, 17), (8, 11), (8, 18), (8, 19), (9, 12), (9, 13), (10, 15), (10, 18), (12, 13), (12, 16), (14, 17), (16, 18)]
            sage: G.show()  # long time

        REFERENCE:

        - [1] Newman, M.E.J. The Structure and function of complex
          networks, SIAM Review vol. 45, no. 2 (2003), pp. 167-256.
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return graph.Graph(networkx.configuration_model([int(i) for i in deg_sequence], seed=seed), loops=True, multiedges=True, sparse=True)

    def DegreeSequenceTree(self, deg_sequence):
        """
        Returns a tree with the given degree sequence. Raises a NetworkX
        error if the proposed degree sequence cannot be that of a tree.

        Since every tree has one more vertex than edge, the degree sequence
        must satisfy len(deg_sequence) - sum(deg_sequence)/2 == 1.

        INPUT:


        -  ``deg_sequence`` - a list of integers with each
           entry corresponding to the expected degree of a different vertex.


        EXAMPLE::

            sage: G = graphs.DegreeSequenceTree([3,1,3,3,1,1,1,2,1])
            sage: G.show()  # long time
        """
        import networkx
        return graph.Graph(networkx.degree_sequence_tree([int(i) for i in deg_sequence]))

    def DegreeSequenceExpected(self, deg_sequence, seed=None):
        """
        Returns a random graph with expected given degree sequence. Raises
        a NetworkX error if the proposed degree sequence cannot be that of
        a graph.

        One requirement is that the sum of the degrees must be even, since
        every edge must be incident with two vertices.

        INPUT:


        -  ``deg_sequence`` - a list of integers with each
           entry corresponding to the expected degree of a different vertex.

        -  ``seed`` - for the random number generator.


        EXAMPLE::

            sage: G = graphs.DegreeSequenceExpected([1,2,3,2,3])
            sage: G.edges(labels=False)
            [(0, 2), (0, 3), (1, 1), (1, 4), (2, 3), (2, 4), (3, 4), (4, 4)]
            sage: G.show()  # long time

        REFERENCE:

        - [1] Chung, Fan and Lu, L. Connected components in random
          graphs with given expected degree
          sequences. Ann. Combinatorics (6), 2002 pp. 125-145.
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return graph.Graph(networkx.expected_degree_graph([int(i) for i in deg_sequence], seed=seed), loops=True)

###########################################################################
#   Graph Iterators
###########################################################################

    def __call__(self, vertices=None, property=lambda x: True, augment='edges',
        size=None, deg_seq=None, degree_sequence=None, loops=False, implementation='c_graph',
        sparse=True, copy = True):
        """
        Accesses the generator of isomorphism class representatives.
        Iterates over distinct, exhaustive representatives. See the docstring
        of this class for full documentation.

        EXAMPLES:

        Print graphs on 3 or less vertices::

            sage: for G in graphs(3, augment='vertices'):
            ...    print G
            Graph on 0 vertices
            Graph on 1 vertex
            Graph on 2 vertices
            Graph on 3 vertices
            Graph on 3 vertices
            Graph on 3 vertices
            Graph on 2 vertices
            Graph on 3 vertices

        ::

            sage: for g in graphs():
            ...    if g.num_verts() > 3: break
            ...    print g
            Graph on 0 vertices
            Graph on 1 vertex
            Graph on 2 vertices
            Graph on 2 vertices
            Graph on 3 vertices
            Graph on 3 vertices
            Graph on 3 vertices
            Graph on 3 vertices

        For more examples, see the class level documentation, or type::

            sage: graphs? # not tested

        REFERENCE:

        - Brendan D. McKay, Isomorph-Free Exhaustive generation.
          Journal of Algorithms Volume 26, Issue 2, February 1998,
          pages 306-324.
        """
        from sage.graphs.all import Graph
        from sage.misc.superseded import deprecation
        from copy import copy as copyfun

        if deg_seq is not None:
            deprecation(11927, "The argument name deg_seq is deprecated. It will be "
                        "removed in a future release of Sage. So, please use "
                        "degree_sequence instead.")
        if degree_sequence is None:
            degree_sequence=deg_seq
        if degree_sequence is not None:
            if vertices is None:
                raise NotImplementedError
            if len(degree_sequence) != vertices or sum(degree_sequence)%2 or sum(degree_sequence) > vertices*(vertices-1):
                raise ValueError("Invalid degree sequence.")
            degree_sequence = sorted(degree_sequence)
            if augment == 'edges':
                property = lambda x: all([degree_sequence[i] >= d for i,d in enumerate(sorted(x.degree()))])
                extra_property = lambda x: degree_sequence == sorted(x.degree())
            else:
                property = lambda x: all([degree_sequence[i] >= d for i,d in enumerate(sorted(x.degree() + [0]*(vertices-x.num_verts()) ))])
                extra_property = lambda x: x.num_verts() == vertices and degree_sequence == sorted(x.degree())
        elif size is not None:
            extra_property = lambda x: x.size() == size
        else:
            extra_property = lambda x: True
        if augment == 'vertices':
            if vertices is None:
                raise NotImplementedError
            g = Graph(loops=loops, implementation=implementation, sparse=sparse)
            for gg in canaug_traverse_vert(g, [], vertices, property, loops=loops, implementation=implementation, sparse=sparse):
                if extra_property(gg):
                    yield copyfun(gg) if copy else gg
        elif augment == 'edges':
            if vertices is None:
                from sage.rings.all import Integer
                vertices = Integer(0)
                while True:
                    for g in self(vertices, loops=loops, implementation=implementation, sparse=sparse):
                        yield copyfun(g) if copy else g
                    vertices += 1
            g = Graph(vertices, loops=loops, implementation=implementation, sparse=sparse)
            gens = []
            for i in range(vertices-1):
                gen = range(i)
                gen.append(i+1); gen.append(i)
                gen += range(i+2, vertices)
                gens.append(gen)
            for gg in canaug_traverse_edge(g, gens, property, loops=loops, implementation=implementation, sparse=sparse):
                if extra_property(gg):
                    yield copyfun(gg) if copy else gg
        else:
            raise NotImplementedError

    def line_graph_forbidden_subgraphs(self):
        r"""
        Returns the 9 forbidden subgraphs of a line graph.

        `Wikipedia article on the line graphs
        <http://en.wikipedia.org/wiki/Line_graph>`_

        The graphs are returned in the ordering given by the Wikipedia
        drawing, read from left to right and from top to bottom.

        EXAMPLE::

            sage: graphs.line_graph_forbidden_subgraphs()
            [Claw graph: Graph on 4 vertices,
            Graph on 6 vertices,
            Graph on 6 vertices,
            Graph on 5 vertices,
            Graph on 6 vertices,
            Graph on 6 vertices,
            Graph on 6 vertices,
            Graph on 6 vertices,
            Graph on 5 vertices]

        """
        from sage.graphs.all import Graph
        graphs = [self.ClawGraph()]

        graphs.append(Graph({
                    0: [1, 2, 3],
                    1: [2, 3],
                    4: [2],
                    5: [3]
                    }))

        graphs.append(Graph({
                    0: [1, 2, 3, 4],
                    1: [2, 3, 4],
                    3: [4],
                    2: [5]
                    }))

        graphs.append(Graph({
                    0: [1, 2, 3],
                    1: [2, 3],
                    4: [2, 3]
                    }))

        graphs.append(Graph({
                    0: [1, 2, 3],
                    1: [2, 3],
                    4: [2],
                    5: [3, 4]
                    }))

        graphs.append(Graph({
                    0: [1, 2, 3, 4],
                    1: [2, 3, 4],
                    3: [4],
                    5: [2, 0, 1]
                    }))

        graphs.append(Graph({
                    5: [0, 1, 2, 3, 4],
                    0: [1, 4],
                    2: [1, 3],
                    3: [4]
                    }))

        graphs.append(Graph({
                    1: [0, 2, 3, 4],
                    3: [0, 4],
                    2: [4, 5],
                    4: [5]
                    }))

        graphs.append(Graph({
                    0: [1, 2, 3],
                    1: [2, 3, 4],
                    2: [3, 4],
                    3: [4]
                    }))

        return graphs

    def trees(self, vertices):
        r"""
        Returns a generator of the distinct trees on a fixed number of vertices.

        INPUT:

        -  ``vertices`` - the size of the trees created.

        OUTPUT:

        A generator which creates an exhaustive, duplicate-free listing
        of the connected free (unlabeled) trees with ``vertices`` number
        of vertices.  A tree is a graph with no cycles.

        ALGORITHM:

        Uses an algorithm that generates each new tree
        in constant time.  See the documentation for, and implementation
        of, the :mod:`sage.graphs.trees` module, including a citation.

        EXAMPLES:

        We create an iterator, then loop over its elements. ::

            sage: tree_iterator = graphs.trees(7)
            sage: for T in tree_iterator:
            ...     print T.degree_sequence()
            [2, 2, 2, 2, 2, 1, 1]
            [3, 2, 2, 2, 1, 1, 1]
            [3, 2, 2, 2, 1, 1, 1]
            [4, 2, 2, 1, 1, 1, 1]
            [3, 3, 2, 1, 1, 1, 1]
            [3, 3, 2, 1, 1, 1, 1]
            [4, 3, 1, 1, 1, 1, 1]
            [3, 2, 2, 2, 1, 1, 1]
            [4, 2, 2, 1, 1, 1, 1]
            [5, 2, 1, 1, 1, 1, 1]
            [6, 1, 1, 1, 1, 1, 1]

        The number of trees on the first few vertex counts.
        This is sequence A000055 in Sloane's OEIS. ::

            sage: [len(list(graphs.trees(i))) for i in range(0, 15)]
            [1, 1, 1, 1, 2, 3, 6, 11, 23, 47, 106, 235, 551, 1301, 3159]
        """
        from trees import TreeIterator
        return iter(TreeIterator(vertices))

    def nauty_geng(self, options="", debug=False):
        r"""
        Returns a generator which creates graphs from nauty's geng program.

        .. note::

            Due to license restrictions, the nauty package is distributed
            as a Sage optional package.  At a system command line, execute
            ``sage -i nauty`` to see the nauty license and install the
            package.

        INPUT:

        - ``options`` - a string passed to  geng  as if it was run at
          a system command line. At a minimum, you *must* pass the
          number of vertices you desire.  Sage expects the graphs to be
          in nauty's "graph6" format, do not set an option to change
          this default or results will be unpredictable.

        - ``debug`` - default: ``False`` - if ``True`` the first line of
          geng's output to standard error is captured and the first call
          to the generator's ``next()`` function will return this line
          as a string.  A line leading with ">A" indicates a successful
          initiation of the program with some information on the arguments,
          while a line beginning with ">E" indicates an error with the input.

        The possible options, obtained as output of ``geng --help``::

                 n    : the number of vertices
            mine:maxe : a range for the number of edges
                        #:0 means '# or more' except in the case 0:0
              res/mod : only generate subset res out of subsets 0..mod-1

                -c    : only write connected graphs
                -C    : only write biconnected graphs
                -t    : only generate triangle-free graphs
                -f    : only generate 4-cycle-free graphs
                -b    : only generate bipartite graphs
                            (-t, -f and -b can be used in any combination)
                -m    : save memory at the expense of time (only makes a
                            difference in the absence of -b, -t, -f and n <= 28).
                -d#   : a lower bound for the minimum degree
                -D#   : a upper bound for the maximum degree
                -v    : display counts by number of edges
                -l    : canonically label output graphs

                -q    : suppress auxiliary output (except from -v)

        Options which cause geng to use an output format different
        than the graph6 format are not listed above (-u, -g, -s, -y, -h)
        as they will confuse the creation of a Sage graph.  The res/mod
        option can be useful when using the output in a routine run
        several times in parallel.

        OUTPUT:

        A generator which will produce the graphs as Sage graphs.
        These will be simple graphs: no loops, no multiple edges, no
        directed edges.

        EXAMPLES:

        The generator can be used to construct graphs for testing,
        one at a time (usually inside a loop).  Or it can be used to
        create an entire list all at once if there is sufficient memory
        to contain it.  ::

            sage: gen = graphs.nauty_geng("2") # optional nauty
            sage: gen.next() # optional nauty
            Graph on 2 vertices
            sage: gen.next() # optional nauty
            Graph on 2 vertices
            sage: gen.next() # optional nauty
            Traceback (most recent call last):
            ...
            StopIteration: Exhausted list of graphs from nauty geng

        A list of all graphs on 7 vertices.  This agrees with
        Sloane's OEIS sequence A000088.  ::

            sage: gen = graphs.nauty_geng("7") # optional nauty
            sage: len(list(gen))  # optional nauty
            1044

        A list of just the connected graphs on 7 vertices.  This agrees with
        Sloane's OEIS sequence A001349.  ::

            sage: gen = graphs.nauty_geng("7 -c") # optional nauty
            sage: len(list(gen))  # optional nauty
            853

        The ``debug`` switch can be used to examine geng's reaction
        to the input in the ``options`` string.  We illustrate success.
        (A failure will be a string beginning with ">E".)  Passing the
        "-q" switch to geng will supress the indicator of a
        successful initiation.  ::

            sage: gen = graphs.nauty_geng("4", debug=True) # optional nauty
            sage: print gen.next() # optional nauty
            >A nauty-geng -d0D3 n=4 e=0-6
        """
        import subprocess
        from sage.misc.package import is_package_installed
        if not is_package_installed("nauty"):
            raise TypeError, "the optional nauty package is not installed"
        sp = subprocess.Popen("nauty-geng {0}".format(options), shell=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, close_fds=True)
        if debug:
            yield sp.stderr.readline()
        gen = sp.stdout
        while True:
            try:
                s = gen.next()
            except StopIteration:
                raise StopIteration("Exhausted list of graphs from nauty geng")
            G = graph.Graph(s[:-1], format='graph6')
            yield G


    def cospectral_graphs(self, vertices, matrix_function=lambda g: g.adjacency_matrix(), graphs=None):
        r"""
        Find all sets of graphs on ``vertices`` vertices (with
        possible restrictions) which are cospectral with respect to a
        constructed matrix.

        INPUT:

        - ``vertices`` - The number of vertices in the graphs to be tested

        - ``matrix_function`` - A function taking a graph and giving back
          a matrix.  This defaults to the adjacency matrix.  The spectra
          examined are the spectra of these matrices.

        - ``graphs`` - One of three things:

           - ``None`` (default) - test all graphs having ``vertices``
             vertices

           - a function taking a graph and returning ``True`` or ``False``
             - test only the graphs on ``vertices`` vertices for which
             the function returns ``True``

           - a list of graphs (or other iterable object) - these graphs
             are tested for cospectral sets.  In this case,
             ``vertices`` is ignored.

        OUTPUT:

           A list of lists of graphs.  Each sublist will be a list of
           cospectral graphs (lists of cadinality 1 being omitted).


        EXAMPLES::

            sage: g=graphs.cospectral_graphs(5)
            sage: sorted(sorted(g.graph6_string() for g in glist) for glist in g)
            [['Dr?', 'Ds_']]
            sage: g[0][1].am().charpoly()==g[0][1].am().charpoly()
            True

        There are two sets of cospectral graphs on six vertices with no isolated vertices::

            sage: g=graphs.cospectral_graphs(6, graphs=lambda x: min(x.degree())>0)
            sage: sorted(sorted(g.graph6_string() for g in glist) for glist in g)
            [['Ep__', 'Er?G'], ['ExGg', 'ExoG']]
            sage: g[0][1].am().charpoly()==g[0][1].am().charpoly()
            True
            sage: g[1][1].am().charpoly()==g[1][1].am().charpoly()
            True

        There is one pair of cospectral trees on eight vertices::

            sage: g=graphs.cospectral_graphs(6, graphs=graphs.trees(8))
            sage: sorted(sorted(g.graph6_string() for g in glist) for glist in g)
            [['GiPC?C', 'GiQCC?']]
            sage: g[0][1].am().charpoly()==g[0][1].am().charpoly()
            True

        There are two sets of cospectral graphs (with respect to the
        Laplacian matrix) on six vertices::

            sage: g=graphs.cospectral_graphs(6, matrix_function=lambda g: g.laplacian_matrix())
            sage: sorted(sorted(g.graph6_string() for g in glist) for glist in g)
            [['Edq_', 'ErcG'], ['Exoo', 'EzcG']]
            sage: g[0][1].laplacian_matrix().charpoly()==g[0][1].laplacian_matrix().charpoly()
            True
            sage: g[1][1].laplacian_matrix().charpoly()==g[1][1].laplacian_matrix().charpoly()
            True

        To find cospectral graphs with respect to the normalized
        Laplacian, assuming the graphs do not have an isolated vertex, it
        is enough to check the spectrum of the matrix `D^{-1}A`, where `D`
        is the diagonal matrix of vertex degrees, and A is the adjacency
        matrix.  We find two such cospectral graphs (for the normalized
        Laplacian) on five vertices::

            sage: def DinverseA(g):
            ...     A=g.adjacency_matrix().change_ring(QQ)
            ...     for i in range(g.order()):
            ...         A.rescale_row(i, 1/len(A.nonzero_positions_in_row(i)))
            ...     return A
            sage: g=graphs.cospectral_graphs(5, matrix_function=DinverseA, graphs=lambda g: min(g.degree())>0)
            sage: sorted(sorted(g.graph6_string() for g in glist) for glist in g)
            [['Dlg', 'Ds_']]
            sage: g[0][1].laplacian_matrix(normalized=True).charpoly()==g[0][1].laplacian_matrix(normalized=True).charpoly()
            True
        """
        from sage.graphs.all import graphs as graph_gen
        if graphs is None:
            graph_list=graph_gen(vertices)
        elif callable(graphs):
            graph_list=iter(g for g in graph_gen(vertices) if graphs(g))
        else:
            graph_list=iter(graphs)

        from collections import defaultdict
        charpolys=defaultdict(list)
        for g in graph_list:
            cp=matrix_function(g).charpoly()
            charpolys[cp].append(g)

        cospectral_graphs=[]
        for cp,g_list in charpolys.items():
            if len(g_list)>1:
                cospectral_graphs.append(g_list)

        return cospectral_graphs


def canaug_traverse_vert(g, aut_gens, max_verts, property, dig=False, loops=False, implementation='c_graph', sparse=True):
    """
    Main function for exhaustive generation. Recursive traversal of a
    canonically generated tree of isomorph free (di)graphs satisfying a
    given property.

    INPUT:


    -  ``g`` - current position on the tree.

    -  ``aut_gens`` - list of generators of Aut(g), in
       list notation.

    -  ``max_verts`` - when to retreat.

    -  ``property`` - check before traversing below g.

    -  ``degree_sequence`` - specify a degree sequence to try to
       obtain.


    EXAMPLES::

        sage: from sage.graphs.graph_generators import canaug_traverse_vert
        sage: list(canaug_traverse_vert(Graph(), [], 3, lambda x: True))
        [Graph on 0 vertices, ... Graph on 3 vertices]

    The best way to access this function is through the graphs()
    iterator:

    Print graphs on 3 or less vertices.

    ::

        sage: for G in graphs(3, augment='vertices'):
        ...    print G
        ...
        Graph on 0 vertices
        Graph on 1 vertex
        Graph on 2 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 2 vertices
        Graph on 3 vertices

    Print digraphs on 2 or less vertices.

    ::

        sage: for D in digraphs(2, augment='vertices'):
        ...    print D
        ...
        Digraph on 0 vertices
        Digraph on 1 vertex
        Digraph on 2 vertices
        Digraph on 2 vertices
        Digraph on 2 vertices
    """
    from sage.graphs.generic_graph_pyx import binary
    from sage.groups.perm_gps.partn_ref.refinement_graphs import search_tree


    if not property(g):
        return
    yield g

    n = g.order()
    if n < max_verts:

        # build a list representing C(g) - the vertex to be added
        # is at the end, so only specify which edges...
        # in the case of graphs, there are n possibilities,
        # and in the case of digraphs, there are 2*n.
        if dig:
            possibilities = 2*n
        else:
            possibilities = n
        num_roots = 2**possibilities
        children = [-1]*num_roots

        # union-find C(g) under Aut(g)
        for gen in aut_gens:
            for i in xrange(len(children)):
                k = 0
                for j in xrange(possibilities):
                    if (1 << j)&i:
                        if dig and j >= n:
                            k += (1 << (gen[j-n]+n))
                        else:
                            k += (1 << gen[j])
                while children[k] != -1:
                    k = children[k]
                while children[i] != -1:
                    i = children[i]
                if i != k:
                    # union i & k
                    smaller, larger = sorted([i,k])
                    children[larger] = smaller
                    num_roots -= 1

        # find representatives of orbits of C(g)
        roots = []
        found_roots = 0
        i = 0
        while found_roots < num_roots:
            if children[i] == -1:
                found_roots += 1
                roots.append(i)
            i += 1
        for i in roots:
            # construct a z for each number in roots...
            z = g.copy(implementation=implementation, sparse=sparse)
            z.add_vertex(n)
            edges = []
            if dig:
                index = 0
                while index < possibilities/2:
                    if (1 << index)&i:
                        edges.append((index,n))
                    index += 1
                while index < possibilities:
                    if (1 << index)&i:
                        edges.append((n,index-n))
                    index += 1
            else:
                index = 0
                while (1 << index) <= i:
                    if (1 << index)&i:
                        edges.append((index,n))
                    index += 1
            z.add_edges(edges)
            z_s = []
            if property(z):
                z_s.append(z)
            if loops:
                z = z.copy(implementation=implementation, sparse=sparse)
                z.add_edge((n,n))
                if property(z):
                    z_s.append(z)
            for z in z_s:
                z_aut_gens, _, canonical_relabeling = search_tree(z, [z.vertices()], certify=True, dig=(dig or loops))
                cut_vert = 0
                while canonical_relabeling[cut_vert] != n:
                    cut_vert += 1
                sub_verts = [v for v in z if v != cut_vert]
                m_z = z.subgraph(sub_verts)

                if m_z == g:
                    for a in canaug_traverse_vert(z, z_aut_gens, max_verts, property, dig=dig, loops=loops, implementation=implementation, sparse=sparse):
                        yield a
                else:
                    for possibility in check_aut(z_aut_gens, cut_vert, n):
                        if m_z.relabel(possibility, inplace=False) == g:
                            for a in canaug_traverse_vert(z, z_aut_gens, max_verts, property, dig=dig, loops=loops, implementation=implementation, sparse=sparse):
                                yield a
                            break

def check_aut(aut_gens, cut_vert, n):
    """
    Helper function for exhaustive generation.

    At the start, check_aut is given a set of generators for the
    automorphism group, aut_gens. We already know we are looking for
    an element of the auto- morphism group that sends cut_vert to n,
    and check_aut generates these for the canaug_traverse function.

    EXAMPLE: Note that the last two entries indicate that none of the
    automorphism group has yet been searched - we are starting at the
    identity [0, 1, 2, 3] and so far that is all we have seen. We
    return automorphisms mapping 2 to 3.

    ::

        sage: from sage.graphs.graph_generators import check_aut
        sage: list( check_aut( [ [0, 3, 2, 1], [1, 0, 3, 2], [2, 1, 0, 3] ], 2, 3))
        [[1, 0, 3, 2], [1, 2, 3, 0]]
    """
    from copy import copy
    perm = range(n+1)
    seen_perms = [perm]
    unchecked_perms = [perm]
    while len(unchecked_perms) != 0:
        perm = unchecked_perms.pop(0)
        for gen in aut_gens:
            new_perm = copy(perm)
            for i in xrange(len(perm)):
                new_perm[i] = gen[perm[i]]
            if new_perm not in seen_perms:
                seen_perms.append(new_perm)
                unchecked_perms.append(new_perm)
                if new_perm[cut_vert] == n:
                    yield new_perm

def canaug_traverse_edge(g, aut_gens, property, dig=False, loops=False, implementation='c_graph', sparse=True):
    """
    Main function for exhaustive generation. Recursive traversal of a
    canonically generated tree of isomorph free graphs satisfying a
    given property.

    INPUT:


    -  ``g`` - current position on the tree.

    -  ``aut_gens`` - list of generators of Aut(g), in
       list notation.

    -  ``property`` - check before traversing below g.


    EXAMPLES::

        sage: from sage.graphs.graph_generators import canaug_traverse_edge
        sage: G = Graph(3)
        sage: list(canaug_traverse_edge(G, [], lambda x: True))
        [Graph on 3 vertices, ... Graph on 3 vertices]

    The best way to access this function is through the graphs()
    iterator:

    Print graphs on 3 or less vertices.

    ::

        sage: for G in graphs(3):
        ...    print G
        ...
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices

    Print digraphs on 3 or less vertices.

    ::

        sage: for G in digraphs(3):
        ...    print G
        ...
        Digraph on 3 vertices
        Digraph on 3 vertices
        ...
        Digraph on 3 vertices
        Digraph on 3 vertices
    """
    from sage.graphs.generic_graph_pyx import binary
    from sage.groups.perm_gps.partn_ref.refinement_graphs import search_tree
    if not property(g):
        return
    yield g
    n = g.order()
    if dig:
        max_size = n*(n-1)
    else:
        max_size = (n*(n-1))>>1 # >> 1 is just / 2 (this is n choose 2)
    if loops: max_size += n
    if g.size() < max_size:
        # build a list representing C(g) - the edge to be added
        # is one of max_size choices
        if dig:
            children = [[(j,i) for i in xrange(n)] for j in xrange(n)]
        else:
            children = [[(j,i) for i in xrange(j)] for j in xrange(n)]
        # union-find C(g) under Aut(g)
        orbits = range(n)
        for gen in aut_gens:
            for iii in xrange(n):
                if orbits[gen[iii]] != orbits[iii]:
                    temp = orbits[gen[iii]]
                    for jjj in xrange(n):
                        if orbits[jjj] == temp:
                            orbits[jjj] = orbits[iii]
                if dig:
                    jjj_range = range(iii) + range(iii+1, n)
                else:
                    jjj_range = xrange(iii) # iii > jjj
                for jjj in jjj_range:
                    i, j = iii, jjj
                    if dig:
                        x, y = gen[i], gen[j]
                    else:
                        y, x = sorted([gen[i], gen[j]])
                    if children[i][j] != children[x][y]:
                        x_val, y_val = x, y
                        i_val, j_val = i, j
                        if dig:
                            while (x_val, y_val) != children[x_val][y_val]:
                                x_val, y_val = children[x_val][y_val]
                            while (i_val, j_val) != children[i_val][j_val]:
                                i_val, j_val = children[i_val][j_val]
                        else:
                            while (x_val, y_val) != children[x_val][y_val]:
                                y_val, x_val = sorted(children[x_val][y_val])
                            while (i_val, j_val) != children[i_val][j_val]:
                                j_val, i_val = sorted(children[i_val][j_val])
                        while (x, y) != (x_val, y_val):
                            xx, yy = x, y
                            x, y = children[x][y]
                            children[xx][yy] = (x_val, y_val)
                        while (i, j) != (i_val, j_val):
                            ii, jj = i, j
                            i, j = children[i][j]
                            children[ii][jj] = (i_val, j_val)
                        if x < i:
                            children[i][j] = (x, y)
                        elif x > i:
                            children[x][y] = (i, j)
                        elif y < j:
                            children[i][j] = (x, y)
                        elif y > j:
                            children[x][y] = (i, j)
                        else:
                            continue
        # find representatives of orbits of C(g)
        roots = []
        for i in range(n):
            if dig:
                j_range = range(i) + range(i+1, n)
            else:
                j_range = range(i)
            for j in j_range:
                if children[i][j] == (i, j):
                    roots.append((i,j))
        if loops:
            seen = []
            for i in xrange(n):
                if orbits[i] not in seen:
                    roots.append((i,i))
                    seen.append(orbits[i])
        for i, j in roots:
            if g.has_edge(i, j):
                continue
            # construct a z for each edge in roots...
            z = g.copy(implementation=implementation, sparse=sparse)
            z.add_edge(i, j)
            if not property(z):
                continue
            z_aut_gens, _, canonical_relabeling = search_tree(z, [z.vertices()], certify=True, dig=(dig or loops))
            relabel_inverse = [0]*n
            for ii in xrange(n):
                relabel_inverse[canonical_relabeling[ii]] = ii
            z_can = z.relabel(canonical_relabeling, inplace=False)
            cut_edge_can = z_can.edges(labels=False, sort=True)[-1]
            cut_edge = [relabel_inverse[cut_edge_can[0]], relabel_inverse[cut_edge_can[1]]]
            if dig:
                cut_edge = tuple(cut_edge)
            else:
                cut_edge = tuple(sorted(cut_edge))

            from copy import copy
            m_z = copy(z)
            m_z.delete_edge(cut_edge)
            if m_z == g:
                for a in canaug_traverse_edge(z, z_aut_gens, property, dig=dig, loops=loops, implementation=implementation, sparse=sparse):
                    yield a
            else:
                for possibility in check_aut_edge(z_aut_gens, cut_edge, i, j, n, dig=dig):
                    if m_z.relabel(possibility, inplace=False) == g:
                        for a in canaug_traverse_edge(z, z_aut_gens, property, dig=dig, loops=loops, implementation=implementation, sparse=sparse):
                            yield a
                        break

def check_aut_edge(aut_gens, cut_edge, i, j, n, dig=False):
    """
    Helper function for exhaustive generation.

    At the start, check_aut_edge is given a set of generators for the
    automorphism group, aut_gens. We already know we are looking for
    an element of the auto- morphism group that sends cut_edge to {i,
    j}, and check_aut generates these for the canaug_traverse
    function.

    EXAMPLE: Note that the last two entries indicate that none of the
    automorphism group has yet been searched - we are starting at the
    identity [0, 1, 2, 3] and so far that is all we have seen. We
    return automorphisms mapping 2 to 3.

    ::

        sage: from sage.graphs.graph_generators import check_aut
        sage: list( check_aut( [ [0, 3, 2, 1], [1, 0, 3, 2], [2, 1, 0, 3] ], 2, 3))
        [[1, 0, 3, 2], [1, 2, 3, 0]]
    """
    from copy import copy
    perm = range(n)
    seen_perms = [perm]
    unchecked_perms = [perm]
    while len(unchecked_perms) != 0:
        perm = unchecked_perms.pop(0)
        for gen in aut_gens:
            new_perm = copy(perm)
            for ii in xrange(n):
                new_perm[ii] = gen[perm[ii]]
            if new_perm not in seen_perms:
                seen_perms.append(new_perm)
                unchecked_perms.append(new_perm)
                if new_perm[cut_edge[0]] == i and new_perm[cut_edge[1]] == j:
                    yield new_perm
                if not dig and new_perm[cut_edge[0]] == j and new_perm[cut_edge[1]] == i:
                    yield new_perm


# Easy access to the graph generators from the command line:
graphs = GraphGenerators()


####################
# Helper functions #
####################

def _circle_embedding(g, vertices, center=(0, 0), radius=1, shift=0):
    r"""
    Set some vertices on a circle in the embedding of a graph G.

    This method modifies the graph's embedding so that the vertices
    listed in ``vertices`` appear in this ordering on a circle of given
    radius and center. The ``shift`` parameter is actually a rotation of
    the circle. A value of ``shift=1`` will replace in the drawing the
    `i`-th element of the list by the `(i-1)`-th. Non-integer values are
    admissible, and a value of `\alpha` corresponds to a rotation of the
    circle by an angle of `\alpha 2\pi/n` (where `n` is the number of
    vertices set on the circle).

    EXAMPLE::

        sage: from sage.graphs.graph_generators import _circle_embedding
        sage: g = graphs.CycleGraph(5)
        sage: _circle_embedding(g, [0, 2, 4, 1, 3], radius=2, shift=.5)
        sage: g.show()
    """
    c_x, c_y = center
    n = len(vertices)
    d = g.get_pos()
    if d is None:
        d = {}

    for i,v in enumerate(vertices):
        i += shift
        v_x = c_x + radius * cos(2*i*pi / n)
        v_y = c_y + radius * sin(2*i*pi / n)
        d[v] = (v_x, v_y)

    g.set_pos(d)

def _line_embedding(g, vertices, first=(0, 0), last=(0, 1)):
    r"""
    Sets some vertices on a line in the embedding of a graph G.

    This method modifies the graph's embedding so that the vertices of
    ``vertices`` appear on a line, where the position of ``vertices[0]``
    is the pair ``first`` and the position of ``vertices[-1]`` is
    ``last``. The vertices are evenly spaced.

    EXAMPLE::

        sage: from sage.graphs.graph_generators import _line_embedding
        sage: g = graphs.PathGraph(5)
        sage: _line_embedding(g, [0, 2, 4, 1, 3], first=(-1, -1), last=(1, 1))
        sage: g.show()
    """
    n = len(vertices) - 1.

    fx, fy = first
    dx = (last[0] - first[0])/n
    dy = (last[1] - first[1])/n

    d = g.get_pos()
    if d is None:
        d = {}

    for v in vertices:
        d[v] = (fx, fy)
        fx += dx
        fy += dy
