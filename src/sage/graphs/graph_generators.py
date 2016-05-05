# -*- coding: utf-8 -*-
r"""
Common Graphs

All graphs in Sage can be built through the ``graphs`` object. In order to
build a complete graph on 15 elements, one can do::

    sage: g = graphs.CompleteGraph(15)

To get a path with 4 vertices, and the house graph::

    sage: p = graphs.PathGraph(4)
    sage: h = graphs.HouseGraph()

More interestingly, one can get the list of all graphs that Sage knows how to
build by typing ``graphs.`` in Sage and then hitting tab.
"""

# This method appends a list of methods to the doc as a 3xN table.

# Here's the point :
#
# we just have to insert the method's name in this file to add it to
# the tab, and in exchange the doc contains a table of width 3 with
# all methods listed, so that the reading order is Column1, then
# Column2, then Column3. Doing this by hand is hell with Sphinx when
# you need to insert a new method inside of the list !

def __append_to_doc(methods):
    global __doc__
    __doc__ += ("\n.. csv-table::\n"
    "    :class: contentstable\n"
    "    :widths: 33, 33, 33\n"
    "    :delim: |\n\n")

    h = (len(methods)+2)//3
    # Reorders the list of methods for horizontal reading, the only one Sphinx understands
    reordered_methods = [0]*3*h
    for i, m in enumerate(methods):
        reordered_methods[3*(i%h)+(i//h)] = m
    methods = reordered_methods

    # Adding the list to the __doc__ string
    wrap_name = lambda x : ":meth:`"+str(x)+" <GraphGenerators."+str(x)+">`" if x else ""
    while methods:
        a = methods.pop(0)
        b = methods.pop(0)
        c = methods.pop(0)
        __doc__ += "    "+wrap_name(a)+" | "+wrap_name(b)+" | "+wrap_name(c)+"\n"

__doc__ += """
**Basic structures**
"""

__append_to_doc(
    ["BullGraph",
     "ButterflyGraph",
     "CircularLadderGraph",
     "ClawGraph",
     "CycleGraph",
     "CompleteBipartiteGraph",
     "CompleteGraph",
     "CompleteMultipartiteGraph",
     "DiamondGraph",
     "EmptyGraph",
     "Grid2dGraph",
     "GridGraph",
     "HouseGraph",
     "HouseXGraph",
     "LadderGraph",
     "LollipopGraph",
     "PathGraph",
     "StarGraph",
     "ToroidalGrid2dGraph",
     "Toroidal6RegularGrid2dGraph"]
    )

__doc__ += """
**Small Graphs**

A small graph is just a single graph and has no parameter influencing
the number of edges or vertices.
"""

__append_to_doc(
    ["Balaban10Cage",
     "Balaban11Cage",
     "BidiakisCube",
     "BiggsSmithGraph",
     "BlanusaFirstSnarkGraph",
     "BlanusaSecondSnarkGraph",
     "BrinkmannGraph",
     "BrouwerHaemersGraph",
     "BuckyBall",
     "CameronGraph",
     "Cell600",
     "Cell120",
     "ChvatalGraph",
     "ClebschGraph",
     "CoxeterGraph",
     "DesarguesGraph",
     "DejterGraph",
     "DoubleStarSnark",
     "DurerGraph",
     "DyckGraph",
     "EllinghamHorton54Graph",
     "EllinghamHorton78Graph",
     "ErreraGraph",
     "F26AGraph",
     "FlowerSnark",
     "FolkmanGraph",
     "FosterGraph",
     "FranklinGraph",
     "FruchtGraph",
     "GoldnerHararyGraph",
     "GossetGraph",
     "GrayGraph",
     "GrotzschGraph",
     "HallJankoGraph",
     "HarborthGraph",
     "HarriesGraph",
     "HarriesWongGraph",
     "HeawoodGraph",
     "HerschelGraph",
     "HigmanSimsGraph",
     "HoffmanGraph",
     "HoffmanSingletonGraph",
     "HoltGraph",
     "HortonGraph",
     "JankoKharaghaniGraph",
     "JankoKharaghaniTonchevGraph",
     "KittellGraph",
     "KrackhardtKiteGraph",
     "Klein3RegularGraph",
     "Klein7RegularGraph",
     "LocalMcLaughlinGraph",
     "LjubljanaGraph",
     "LivingstoneGraph",
     "M22Graph",
     "MarkstroemGraph",
     "MathonStronglyRegularGraph",
     "McGeeGraph",
     "McLaughlinGraph",
     "MeredithGraph",
     "MoebiusKantorGraph",
     "MoserSpindle",
     "NauruGraph",
     "PappusGraph",
     "PoussinGraph",
     "PerkelGraph",
     "PetersenGraph",
     "RobertsonGraph",
     "SchlaefliGraph",
     "ShrikhandeGraph",
     "SimsGewirtzGraph",
     "SousselierGraph",
     "SylvesterGraph",
     "SzekeresSnarkGraph",
     "ThomsenGraph",
     "TietzeGraph",
     "TruncatedIcosidodecahedralGraph",
     "TruncatedTetrahedralGraph",
     "Tutte12Cage",
     "TutteCoxeterGraph",
     "TutteGraph",
     "WagnerGraph",
     "WatkinsSnarkGraph",
     "WellsGraph",
     "WienerArayaGraph",
     "SuzukiGraph"])

__doc__ += """
**Platonic solids** (ordered ascending by number of vertices)
"""

__append_to_doc(
    ["TetrahedralGraph",
     "OctahedralGraph",
     "HexahedralGraph",
     "IcosahedralGraph",
     "DodecahedralGraph"])

__doc__ += """
**Families of graphs**

A family of graph is an infinite set of graphs which can be indexed by fixed
number of parameters, e.g. two integer parameters. (A method whose name starts
with a small letter does not return a single graph object but a graph iterator
or a list of graphs or ...)
"""

__append_to_doc(
    ["BalancedTree",
     "BarbellGraph",
     "BubbleSortGraph",
     "chang_graphs",
     "CirculantGraph",
     "cospectral_graphs",
     "CubeGraph",
     "DorogovtsevGoltsevMendesGraph",
     "FibonacciTree",
     "FoldedCubeGraph",
     "FriendshipGraph",
     "fullerenes",
     "fusenes",
     "FuzzyBallGraph",
     "GeneralizedPetersenGraph",
     "GoethalsSeidelGraph",
     "HanoiTowerGraph",
     "HararyGraph",
     "HyperStarGraph",
     "JohnsonGraph",
     "KneserGraph",
     "LCFGraph",
     "line_graph_forbidden_subgraphs",
     "MathonPseudocyclicMergingGraph",
     "MathonPseudocyclicStronglyRegularGraph",
     "MycielskiGraph",
     "MycielskiStep",
     "NKStarGraph",
     "NStarGraph",
     "OddGraph",
     "PaleyGraph",
     "PasechnikGraph",
     "petersen_family",
     "planar_graphs",
     "quadrangulations",
     "RingedTree",
     "SierpinskiGasketGraph",
     "SquaredSkewHadamardMatrixGraph",
     "SwitchedSquaredSkewHadamardMatrixGraph",
     "strongly_regular_graph",
     "trees",
     "triangulations",
     "TuranGraph",
     "WheelGraph"])


__doc__ += """
**Graphs from classical geometries over finite fields**

A number of classes of graphs related to geometries over finite fields and
quadrics and Hermitean varieties there.
"""

__append_to_doc(
    ["AffineOrthogonalPolarGraph",
     "AhrensSzekeresGeneralizedQuadrangleGraph",
     "NonisotropicOrthogonalPolarGraph",
     "NonisotropicUnitaryPolarGraph",
     "OrthogonalPolarGraph",
     "SymplecticDualPolarGraph",
     "SymplecticPolarGraph",
     "TaylorTwographDescendantSRG",
     "TaylorTwographSRG",
     "T2starGeneralizedQuadrangleGraph",
     "HaemersGraph",
     "CossidentePenttilaGraph",
     "UnitaryDualPolarGraph",
     "UnitaryPolarGraph"])

__doc__ += """
**Chessboard Graphs**
"""

__append_to_doc(
    ["BishopGraph",
     "KingGraph",
     "KnightGraph",
     "QueenGraph",
     "RookGraph"])

__doc__ += """
**Intersection graphs**

These graphs are generated by geometric representations. The objects of
the representation correspond to the graph vertices and the intersections
of objects yield the graph edges.
"""

__append_to_doc(
    ["IntersectionGraph",
     "IntervalGraph",
     "OrthogonalArrayBlockGraph",
     "PermutationGraph",
     "ToleranceGraph"])

__doc__ += """
**Random graphs**
"""

__append_to_doc(
    ["RandomBarabasiAlbert",
     "RandomBipartite",
     "RandomBoundedToleranceGraph",
     "RandomGNM",
     "RandomGNP",
     "RandomHolmeKim",
     "RandomIntervalGraph",
     "RandomLobster",
     "RandomNewmanWattsStrogatz",
     "RandomRegular",
     "RandomShell",
     "RandomToleranceGraph",
     "RandomTree",
     "RandomTreePowerlaw",
     "RandomTriangulation"])

__doc__ += """
**Graphs with a given degree sequence**
"""

__append_to_doc(
    ["DegreeSequence",
     "DegreeSequenceBipartite",
     "DegreeSequenceConfigurationModel",
     "DegreeSequenceExpected",
     "DegreeSequenceTree"])

__doc__ += """
**Miscellaneous**
"""

__append_to_doc(
    ["WorldMap"]
    )

__doc__ += """

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

- Ed Scheinerman (2010-11-15): added RandomTriangulation

- Minh Van Nguyen (2010-11-26): added more named graphs

- Keshav Kini (2011-02-16): added Shrikhande and Dyck graphs

- David Coudert (2012-02-10): new RandomGNP generator

- David Coudert (2012-08-02): added chessboard graphs: Queen, King,
  Knight, Bishop, and Rook graphs

- Nico Van Cleemput (2013-05-26): added fullerenes

- Nico Van Cleemput (2013-07-01): added benzenoids

- Birk Eisermann (2013-07-29): new section 'intersection graphs',
  added (random, bounded) tolerance graphs

- Marco Cognetta (2016-03-03): added TuranGraph


Functions and methods
---------------------
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

# import from Sage library
import graph
import sage.graphs.strongly_regular_db

class GraphGenerators():
    r"""
    A class consisting of constructors for several common graphs, as well as
    orderly generation of isomorphism class representatives. See the
    :mod:`module's help <sage.graphs.graph_generators>` for a list of supported
    constructors.

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

    Also: see the use of the nauty package for generating graphs
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
    :oeis:`A033995`)

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

    Generate graphs on the fly: (see :oeis:`A000088`)

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

    Generate all simple graphs, allowing loops: (see :oeis:`A000666`)

    ::

        sage: L = list(graphs(5,augment='vertices',loops=True))               # long time
        sage: for i in [0..5]: print i, len([g for g in L if g.order() == i]) # long time
        0 1
        1 2
        2 6
        3 20
        4 90
        5 544

    Generate all graphs with a specified degree sequence (see :oeis:`A002851`)::

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
    survives repeated vertex removal (:trac:`8458`)::

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

###########################################################################
#   Graph Iterators
###########################################################################

    def __call__(self, vertices=None, property=lambda x: True, augment='edges',
        size=None, degree_sequence=None, loops=False, implementation='c_graph',
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


    def nauty_geng(self, options="", debug=False):
        r"""
        Returns a generator which creates graphs from nauty's geng program.

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

        .. SEEALSO::

            :meth:`Graph.is_strongly_regular` -- tests whether a graph is
            strongly regular and/or returns its parameters.

        EXAMPLES:

        The generator can be used to construct graphs for testing,
        one at a time (usually inside a loop).  Or it can be used to
        create an entire list all at once if there is sufficient memory
        to contain it.  ::

            sage: gen = graphs.nauty_geng("2")
            sage: next(gen)
            Graph on 2 vertices
            sage: next(gen)
            Graph on 2 vertices
            sage: next(gen)
            Traceback (most recent call last):
            ...
            StopIteration: Exhausted list of graphs from nauty geng

        A list of all graphs on 7 vertices.  This agrees with
        :oeis:`A000088`.  ::

            sage: gen = graphs.nauty_geng("7")
            sage: len(list(gen))
            1044

        A list of just the connected graphs on 7 vertices.  This agrees with
        :oeis:`A001349`.  ::

            sage: gen = graphs.nauty_geng("7 -c")
            sage: len(list(gen))
            853

        The ``debug`` switch can be used to examine geng's reaction
        to the input in the ``options`` string.  We illustrate success.
        (A failure will be a string beginning with ">E".)  Passing the
        "-q" switch to geng will supress the indicator of a
        successful initiation.  ::

            sage: gen = graphs.nauty_geng("4", debug=True)
            sage: print next(gen)
            >A geng -d0D3 n=4 e=0-6
        """
        import subprocess
        sp = subprocess.Popen("geng {0}".format(options), shell=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, close_fds=True)
        if debug:
            yield sp.stderr.readline()
        gen = sp.stdout
        while True:
            try:
                s = next(gen)
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


        .. SEEALSO::

            :meth:`Graph.is_strongly_regular` -- tests whether a graph is
            strongly regular and/or returns its parameters.

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

    def _read_planar_code(self, code_input):
        r"""
        Returns a generator for the plane graphs in planar code format in
        the file code_input (see [plantri-guide]_).

        A file with planar code starts with a header ``>>planar_code<<``.
        After the header each graph is stored in the following way :

        The first character is the number of vertices, followed by
        n11,...,n1k,null character,n21,...,n2k',null character, ...

        where the n1* are all neighbors of n1 and all n2* are the 
        neighbors of n2, ...
        Besides, these neighbors are enumerated in clockwise order.

        INPUT:

        - ``code_input`` - a file containing valid planar code data.

        OUTPUT:

        A generator which will produce the plane graphs as Sage graphs
        with an embedding set. These will be simple graphs: no loops, no
        multiple edges, no directed edges (unless plantri is asked to give
        the dual graphs instead).

        .. SEEALSO::

            - :meth:`~sage.graphs.generic_graph.GenericGraph.set_embedding`,
              :meth:`~sage.graphs.generic_graph.GenericGraph.get_embedding` --
              get/set methods for embeddings.

        EXAMPLES:

        The following example creates a small planar code file in memory and
        reads it using the ``_read_planar_code`` method:  ::

            sage: import StringIO
            sage: code_input = StringIO.StringIO('>>planar_code<<')
            sage: code_input.write('>>planar_code<<')
            sage: for c in [4,2,3,4,0,1,4,3,0,1,2,4,0,1,3,2,0]:
            ....:     code_input.write('{:c}'.format(c))
            sage: code_input.seek(0)
            sage: gen = graphs._read_planar_code(code_input)
            sage: l = list(gen)
            sage: l
            [Graph on 4 vertices]
            sage: l[0].is_isomorphic(graphs.CompleteGraph(4))
            True
            sage: l[0].get_embedding()
            {1: [2, 3, 4],
             2: [1, 4, 3],
             3: [1, 2, 4],
             4: [1, 3, 2]}

        REFERENCE:

        .. [plantri-guide] http://cs.anu.edu.au/~bdm/plantri/plantri-guide.txt
        """
        #start of code to read planar code

        header = code_input.read(15)
        assert header == '>>planar_code<<', 'Not a valid planar code header'

        #read graph per graph
        while True:
            c = code_input.read(1)
            if len(c)==0:
                return

            # Each graph is stored in the following way :
            #
            # The first character is the number of vertices, followed by
            # n11,...,n1k,null character,n21,...,n2k',null character, ...
            #
            # where the n1* are all neighbors of n1 and all n2* are the
            # neighbors of n2, ...
            #
            # Besides, these neighbors are enumerated in clockwise order.
            order = ord(c)

            zeroCount = 0

            g = [[] for i in range(order)]

            while zeroCount < order:
                c = code_input.read(1)
                if ord(c) == 0:
                    zeroCount += 1
                else:
                    g[zeroCount].append(ord(c))

            # construct graph based on g

            # first taking care that every edge is given twice
            edges_g = {i + 1: [j for j in di if j < i + 1]
                       for i, di in enumerate(g)}

            # then adding half of the loops (if any)
            has_loops = False
            for i, di in enumerate(g):
                Ni = di.count(i + 1)
                if Ni > 1:
                    edges_g[i + 1] += [i + 1] * (Ni / 2)
                    has_loops = True
            G = graph.Graph(edges_g, loops=has_loops)

            if not(G.has_multiple_edges() or has_loops):
                embed_g = {i + 1: di for i, di in enumerate(g)}
                G.set_embedding(embed_g)
            yield(G)

    def fullerenes(self, order, ipr=False):
        r"""
        Returns a generator which creates fullerene graphs using
        the buckygen generator (see [buckygen]_).

        INPUT:

        - ``order`` - a positive even integer smaller than or equal to 254.
          This specifies the number of vertices in the generated fullerenes.

        - ``ipr`` - default: ``False`` - if ``True`` only fullerenes that
          satisfy the Isolated Pentagon Rule are generated. This means that
          no pentagonal faces share an edge.

        OUTPUT:

        A generator which will produce the fullerene graphs as Sage graphs
        with an embedding set. These will be simple graphs: no loops, no
        multiple edges, no directed edges.

        .. SEEALSO::

            - :meth:`~sage.graphs.generic_graph.GenericGraph.set_embedding`,
              :meth:`~sage.graphs.generic_graph.GenericGraph.get_embedding` --
              get/set methods for embeddings.

        EXAMPLES:

        There are 1812 isomers of `\textrm{C}_{60}`, i.e., 1812 fullerene graphs
        on 60 vertices:  ::

            sage: gen = graphs.fullerenes(60)  # optional buckygen
            sage: len(list(gen))  # optional buckygen
            1812

        However, there is only one IPR fullerene graph on 60 vertices: the famous
        Buckminster Fullerene:  ::

            sage: gen = graphs.fullerenes(60, ipr=True)  # optional buckygen
            sage: next(gen)  # optional buckygen
            Graph on 60 vertices
            sage: next(gen)  # optional buckygen
            Traceback (most recent call last):
            ...
            StopIteration

        The unique fullerene graph on 20 vertices is isomorphic to the dodecahedron
        graph. ::

            sage: gen = graphs.fullerenes(20)  # optional buckygen
            sage: g = next(gen)  # optional buckygen
            sage: g.is_isomorphic(graphs.DodecahedralGraph()) # optional buckygen
            True
            sage: g.get_embedding()  # optional buckygen
            {1: [2, 3, 4],
             2: [1, 5, 6],
             3: [1, 7, 8],
             4: [1, 9, 10],
             5: [2, 10, 11],
             6: [2, 12, 7],
             7: [3, 6, 13],
             8: [3, 14, 9],
             9: [4, 8, 15],
             10: [4, 16, 5],
             11: [5, 17, 12],
             12: [6, 11, 18],
             13: [7, 18, 14],
             14: [8, 13, 19],
             15: [9, 19, 16],
             16: [10, 15, 17],
             17: [11, 16, 20],
             18: [12, 20, 13],
             19: [14, 20, 15],
             20: [17, 19, 18]}
            sage: g.plot3d(layout='spring')  # optional buckygen
            Graphics3d Object

        REFERENCE:

        .. [buckygen] \G. Brinkmann, J. Goedgebeur and B.D. McKay, Generation of Fullerenes,
          Journal of Chemical Information and Modeling, 52(11):2910-2918, 2012.
        """
        from sage.misc.package import is_package_installed
        if not is_package_installed("buckygen"):
            raise TypeError("the optional buckygen package is not installed")

        # number of vertices should be positive
        if order < 0:
            raise ValueError("Number of vertices should be positive.")

        # buckygen can only output fullerenes on up to 254 vertices
        if order > 254:
            raise ValueError("Number of vertices should be at most 254.")

        # fullerenes only exist for an even number of vertices, larger than 20
        # and different from 22
        if order % 2 == 1 or order < 20 or order == 22:
            return

        command = 'buckygen -'+('I' if ipr else '')+'d {0}d'.format(order)

        import subprocess
        sp = subprocess.Popen(command, shell=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, close_fds=True)

        for G in graphs._read_planar_code(sp.stdout):
            yield(G)

    def fusenes(self, hexagon_count, benzenoids=False):
        r"""
        Returns a generator which creates fusenes and benzenoids using
        the benzene generator (see [benzene]_). Fusenes are planar
        polycyclic hydrocarbons with all bounded faces hexagons. Benzenoids
        are fusenes that are subgraphs of the hexagonal lattice.

        INPUT:

        - ``hexagon_count`` - a positive integer smaller than or equal to 30.
          This specifies the number of hexagons in the generated benzenoids.

        - ``benzenoids`` - default: ``False`` - if ``True`` only benzenoids are
          generated.

        OUTPUT:

        A generator which will produce the fusenes as Sage graphs
        with an embedding set. These will be simple graphs: no loops, no
        multiple edges, no directed edges.

        .. SEEALSO::

            - :meth:`~sage.graphs.generic_graph.GenericGraph.set_embedding`,
              :meth:`~sage.graphs.generic_graph.GenericGraph.get_embedding` --
              get/set methods for embeddings.

        EXAMPLES:

        There is a unique fusene with 2 hexagons:  ::

            sage: gen = graphs.fusenes(2)  # optional benzene
            sage: len(list(gen))  # optional benzene
            1

        This fusene is naphtalene (`\textrm{C}_{10}\textrm{H}_{8}`).
        In the fusene graph the H-atoms are not stored, so this is
        a graph on just 10 vertices:  ::

            sage: gen = graphs.fusenes(2)  # optional benzene
            sage: next(gen)  # optional benzene
            Graph on 10 vertices
            sage: next(gen)  # optional benzene
            Traceback (most recent call last):
            ...
            StopIteration

        There are 6505 benzenoids with 9 hexagons:  ::

            sage: gen = graphs.fusenes(9, benzenoids=True)  # optional benzene
            sage: len(list(gen))  # optional benzene
            6505

        REFERENCE:

        .. [benzene] \G. Brinkmann, G. Caporossi and P. Hansen, A Constructive Enumeration of Fusenes and Benzenoids,
          Journal of Algorithms, 45:155-166, 2002.
        """
        from sage.misc.package import is_package_installed
        if not is_package_installed("benzene"):
            raise TypeError("the optional benzene package is not installed")

        # number of hexagons should be positive
        if hexagon_count < 0:
            raise ValueError("Number of hexagons should be positive.")

        # benzene is only built for fusenes with up to 30 hexagons
        if hexagon_count > 30:
            raise ValueError("Number of hexagons should be at most 30.")

        # there are no fusenes with 0 hexagons
        if hexagon_count == 0:
            return

        # there is only one unique fusene with 1 hexagon (and benzene doesn't generate it)
        if hexagon_count == 1:
            g = {1:[6, 2], 2:[1, 3], 3:[2, 4], 4:[3, 5], 5:[4, 6], 6:[5, 1]}
            G = graph.Graph(g)
            G.set_embedding(g)
            yield(G)
            return

        command = 'benzene '+('b' if benzenoids else '')+' {0} p'.format(hexagon_count)

        import subprocess
        sp = subprocess.Popen(command, shell=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, close_fds=True)

        for G in graphs._read_planar_code(sp.stdout):
            yield(G)

    def planar_graphs(self, order, minimum_degree=None,
                      minimum_connectivity=None,
                      exact_connectivity=False, only_bipartite=False,
                      dual=False):
        r"""
        An iterator over connected planar graphs using the plantri generator.

        This uses the plantri generator (see [plantri]_) which is available
        through the optional package plantri.

        .. NOTE::

            The non-3-connected graphs will be returned several times, with all
            its possible embeddings.

        INPUT:

        - ``order`` - a positive integer smaller than or equal to 64.
          This specifies the number of vertices in the generated graphs.

        - ``minimum_degree`` - default: ``None`` - a value `\geq 1` and `\leq
          5`, or ``None``. This specifies the minimum degree of the generated
          graphs. If this is ``None`` and the order is 1, then this is set to
          0. If this is ``None`` and the minimum connectivity is specified, then
          this is set to the same value as the minimum connectivity.  If the
          minimum connectivity is also equal to ``None``, then this is set to 1.

        - ``minimum_connectivity`` - default: ``None`` - a value `\geq 1`
          and `\leq 3`, or ``None``. This specifies the minimum connectivity of the
          generated graphs. If this is ``None`` and the minimum degree is
          specified, then this is set to the minimum of the minimum degree
          and 3. If the minimum degree is also equal to ``None``, then this
          is set to 1.

        - ``exact_connectivity`` - default: ``False`` - if ``True`` only
          graphs with exactly the specified connectivity will be generated.
          This option cannot be used with ``minimum_connectivity=3``, or if
          the minimum connectivity is not explicitely set.

        - ``only_bipartite`` - default: ``False`` - if ``True`` only bipartite
          graphs will be generated. This option cannot be used for graphs with
          a minimum degree larger than 3.

        - ``dual`` - default: ``False`` - if ``True`` return instead the
          planar duals of the generated graphs.

        OUTPUT:

        An iterator which will produce all planar graphs with the given
        number of vertices as Sage graphs with an embedding set. These will be
        simple graphs (no loops, no multiple edges, no directed edges)
        unless the option ``dual=True`` is used.

        .. SEEALSO::

            - :meth:`~sage.graphs.generic_graph.GenericGraph.set_embedding`,
              :meth:`~sage.graphs.generic_graph.GenericGraph.get_embedding` --
              get/set methods for embeddings.

        EXAMPLES:

        There are 6 planar graphs on 4 vertices::

            sage: gen = graphs.planar_graphs(4)  # optional plantri
            sage: len(list(gen))  # optional plantri
            6

        Three of these planar graphs are bipartite::

            sage: gen = graphs.planar_graphs(4, only_bipartite=True)  # optional plantri
            sage: len(list(gen))  # optional plantri
            3

        Setting ``dual=True`` gives the planar dual graphs::

            sage: gen = graphs.planar_graphs(4, dual=True)  # optional plantri
            sage: [u for u in list(gen)]  # optional plantri
            [Graph on 4 vertices,
            Multi-graph on 3 vertices,
            Multi-graph on 2 vertices,
            Looped multi-graph on 2 vertices,
            Looped multi-graph on 1 vertex,
            Looped multi-graph on 1 vertex]

        The cycle of length 4 is the only 2-connected bipartite planar graph
        on 4 vertices::

            sage: l = list(graphs.planar_graphs(4, minimum_connectivity=2, only_bipartite=True))  # optional plantri
            sage: l[0].get_embedding()  # optional plantri
            {1: [2, 3],
             2: [1, 4],
             3: [1, 4],
             4: [2, 3]}

        There is one planar graph with one vertex. This graph obviously has
        minimum degree equal to 0::

            sage: list(graphs.planar_graphs(1))  # optional plantri
            [Graph on 1 vertex]
            sage: list(graphs.planar_graphs(1, minimum_degree=1))  # optional plantri
            []

        TESTS:

        The number of edges in a planar graph is equal to the number of edges in
        its dual::

            sage: planar      = list(graphs.planar_graphs(5,dual=True))  # optional -- plantri
            sage: dual_planar = list(graphs.planar_graphs(5,dual=False)) # optional -- plantri
            sage: planar_sizes      = [g.size() for g in planar]         # optional -- plantri
            sage: dual_planar_sizes = [g.size() for g in dual_planar]    # optional -- plantri
            sage: planar_sizes == dual_planar_sizes                      # optional -- plantri
            True

        REFERENCE:

        .. [plantri] \G. Brinkmann and B.D. McKay, Fast generation of planar graphs,
           MATCH-Communications in Mathematical and in Computer Chemistry, 58(2):323-357, 2007.
        """
        from sage.misc.package import is_package_installed
        if not is_package_installed("plantri"):
            raise TypeError("the optional plantri package is not installed")

        # number of vertices should be positive
        if order < 0:
            raise ValueError("Number of vertices should be positive.")

        # plantri can only output general planar graphs on up to 64 vertices
        if order > 64:
            raise ValueError("Number of vertices should be at most 64.")

        if exact_connectivity and minimum_connectivity is None:
            raise ValueError("Minimum connectivity must be specified to use the exact_connectivity option.")

        # minimum connectivity should be None or a number between 1 and 3
        if minimum_connectivity is  not None and not (1 <= minimum_connectivity <= 3):
            raise ValueError("Minimum connectivity should be a number between 1 and 3.")

        # minimum degree should be None or a number between 1 and 5
        if minimum_degree == 0:
            if order != 1:
                raise ValueError("Minimum degree equal to 0 is only possible if the graphs have 1 vertex.")
        elif minimum_degree is not None and not (1 <= minimum_degree <= 5):
            raise ValueError("Minimum degree should be a number between 1 and 5 if the order is greater than 1.")
        elif minimum_degree is None and order == 1:
            minimum_degree = 0

        # check combination of values of minimum degree and minimum connectivity
        if minimum_connectivity is None:
            if minimum_degree is not None:
                minimum_connectivity = min(3, minimum_degree)
            elif minimum_degree is None:
                minimum_degree, minimum_connectivity = 1, 1
        else:
            if minimum_degree is None:
                minimum_degree = minimum_connectivity
            elif (minimum_degree < minimum_connectivity and
                  minimum_degree > 0):
                raise ValueError("Minimum connectivity can be at most the minimum degree.")

        #exact connectivity is not implemented for minimum connectivity 3
        if exact_connectivity and minimum_connectivity==3:
            raise NotImplementedError("Generation of planar graphs with connectivity exactly 3 is not implemented.")

        if only_bipartite and minimum_degree > 3:
            raise NotImplementedError("Generation of bipartite planar graphs with minimum degree 4 or 5 is not implemented.")

        if order == 0:
            return

        minimum_order = {0:1, 1:2, 2:3, 3:4, 4:6, 5:12}[minimum_degree]

        if order < minimum_order:
            return

        if order == 1:
            if minimum_degree == 0:
                G = graph.Graph(1)
                G.set_embedding({0: []})
                yield(G)
            return

        cmd = 'plantri -p{}m{}c{}{}{} {}'
        command = cmd.format('b' if only_bipartite else '',
                             minimum_degree,
                             minimum_connectivity,
                             'x' if exact_connectivity else '',
                             'd' if dual else '',
                             order)

        import subprocess
        sp = subprocess.Popen(command, shell=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, close_fds=True)

        for G in graphs._read_planar_code(sp.stdout):
            yield(G)

    def triangulations(self, order, minimum_degree=None, minimum_connectivity=None,
                       exact_connectivity=False, only_eulerian=False, dual=False):
        r"""
        An iterator over connected planar triangulations using the plantri generator.

        This uses the plantri generator (see [plantri]_) which is available
        through the optional package plantri.

        INPUT:

        - ``order`` - a positive integer smaller than or equal to 64.
          This specifies the number of vertices in the generated triangulations.

        - ``minimum_degree`` - default: ``None`` - a value `\geq 3` and `\leq 5`,
          or ``None``. This specifies the minimum degree of the generated
          triangulations. If this is ``None`` and the minimum connectivity
          is specified, then this is set to the same value as the minimum
          connectivity. If the minimum connectivity is also equal to ``None``,
          then this is set to 3.

        - ``minimum_connectivity`` - default: ``None`` - a value `\geq 3` and
          `\leq 5`, or ``None``. This specifies the minimum connectivity of the
          generated triangulations. If this is ``None`` and the minimum degree
          is specified, then this is set to the minimum of the minimum degree
          and 3. If the minimum degree is also equal to ``None``, then this is
          set to 3.

        - ``exact_connectivity`` - default: ``False`` - if ``True`` only
          triangulations with exactly the specified connectivity will be generated.
          This option cannot be used with ``minimum_connectivity=3``, or if
          the minimum connectivity is not explicitely set.

        - ``only_eulerian`` - default: ``False`` - if ``True`` only Eulerian
          triangulations will be generated. This option cannot be used if the
          minimum degree is explicitely set to anything else than 4.

        - ``dual`` - default: ``False`` - if ``True`` return instead the
          planar duals of the generated graphs.

        OUTPUT:

        An iterator which will produce all planar triangulations with the given
        number of vertices as Sage graphs with an embedding set. These will be
        simple graphs (no loops, no multiple edges, no directed edges).

        .. SEEALSO::

            - :meth:`~sage.graphs.generic_graph.GenericGraph.set_embedding`,
              :meth:`~sage.graphs.generic_graph.GenericGraph.get_embedding` --
              get/set methods for embeddings.

            - :meth:`~sage.graphs.graph_generators.GraphGenerators.RandomTriangulation`
              -- build a random triangulation.

        EXAMPLES:

        The unique planar embedding of the `K_4` is the only planar triangulations
        on 4 vertices::

            sage: gen = graphs.triangulations(4)    # optional plantri
            sage: [g.get_embedding() for g in gen]  # optional plantri
            [{1: [2, 3, 4], 2: [1, 4, 3], 3: [1, 2, 4], 4: [1, 3, 2]}]

        but, of course, this graph is not Eulerian::

            sage: gen = graphs.triangulations(4, only_eulerian=True)  # optional plantri
            sage: len(list(gen))                                      # optional plantri
            0

        The unique Eulerian triangulation on 6 vertices is isomorphic to the octahedral
        graph. ::

            sage: gen = graphs.triangulations(6, only_eulerian=True)  # optional plantri
            sage: g = next(gen)                                       # optional plantri
            sage: g.is_isomorphic(graphs.OctahedralGraph())           # optional plantri
            True

        An overview of the number of 5-connected triangulations on up to 22 vertices. This
        agrees with :oeis:`A081621`::

            sage: for i in range(12, 23):                                             # optional plantri
            ....:     L = len(list(graphs.triangulations(i, minimum_connectivity=5))) # optional plantri
            ....:     print("{}   {:3d}".format(i,L))                                 # optional plantri
            12     1
            13     0
            14     1
            15     1
            16     3
            17     4
            18    12
            19    23
            20    71
            21   187
            22   627

        The minimum connectivity can be at most the minimum degree::

            sage: gen = next(graphs.triangulations(10, minimum_degree=3, minimum_connectivity=5))  # optional plantri
            Traceback (most recent call last):
            ...
            ValueError: Minimum connectivity can be at most the minimum degree.

        There are 5 triangulations with 9 vertices and minimum degree equal to 4
        that are 3-connected, but only one of them is not 4-connected::

            sage: len([g for g in graphs.triangulations(9, minimum_degree=4, minimum_connectivity=3)]) # optional plantri
            5
            sage: len([g for g in graphs.triangulations(9, minimum_degree=4, minimum_connectivity=3, exact_connectivity=True)]) # optional plantri
            1

        Setting ``dual=True`` gives the planar dual graphs::

            sage: [len(g) for g in graphs.triangulations(9, minimum_degree=4, minimum_connectivity=3, dual=True)]  # optional plantri
            [14, 14, 14, 14, 14]

        TESTS::

            sage: [g.size() for g in graphs.triangulations(6, minimum_connectivity=3)] # optional plantri
            [12, 12]
        """
        from sage.misc.package import is_package_installed
        if not is_package_installed("plantri"):
            raise TypeError("the optional plantri package is not installed")

        # number of vertices should be positive
        if order < 0:
            raise ValueError("Number of vertices should be positive.")

        # plantri can only output planar triangulations on up to 64 vertices
        if order > 64:
            raise ValueError("Number of vertices should be at most 64.")

        if exact_connectivity and minimum_connectivity is None:
            raise ValueError("Minimum connectivity must be specified to use the exact_connectivity option.")

        # minimum connectivity should be None or a number between 3 and 5
        if minimum_connectivity is  not None and not (3 <= minimum_connectivity <= 5):
            raise ValueError("Minimum connectivity should be None or a number between 3 and 5.")

        # minimum degree should be None or a number between 3 and 5
        if minimum_degree is  not None and not (3 <= minimum_degree <= 5):
            raise ValueError("Minimum degree should be None or a number between 3 and 5.")

        # for Eulerian triangulations the minimum degree is set to 4 (unless it was already specifically set)
        if only_eulerian and minimum_degree is None:
            minimum_degree = 4

        # check combination of values of minimum degree and minimum connectivity
        if minimum_connectivity is None:
            if minimum_degree is not None:
                minimum_connectivity = min(3, minimum_degree)
            else:
                minimum_degree, minimum_connectivity = 3, 3
        else:
            if minimum_degree is None:
                minimum_degree = minimum_connectivity
            elif minimum_degree < minimum_connectivity:
                raise ValueError("Minimum connectivity can be at most the minimum degree.")

        #exact connectivity is not implemented for minimum connectivity equal to minimum degree
        if exact_connectivity and minimum_connectivity==minimum_degree:
            raise NotImplementedError("Generation of triangulations with minimum connectivity equal to minimum degree is not implemented.")

        minimum_order = {3:4, 4:6, 5:12}[minimum_degree]

        if order < minimum_order:
            return

        if only_eulerian and order < 6:
            return

        cmd = 'plantri -{}m{}c{}{}{} {}'
        command = cmd.format('b' if only_eulerian else '',
                             minimum_degree,
                             minimum_connectivity,
                             'x' if exact_connectivity else '',
                             'd' if dual else '',
                             order)

        import subprocess
        sp = subprocess.Popen(command, shell=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, close_fds=True)

        for G in graphs._read_planar_code(sp.stdout):
            yield(G)

    def quadrangulations(self, order, minimum_degree=None, minimum_connectivity=None,
                         no_nonfacial_quadrangles=False, dual=False):
        r"""
        An iterator over planar quadrangulations using the plantri generator.

        This uses the plantri generator (see [plantri]_) which is available
        through the optional package plantri.

        INPUT:

        - ``order`` - a positive integer smaller than or equal to 64.
          This specifies the number of vertices in the generated quadrangulations.

        - ``minimum_degree`` - default: ``None`` - a value `\geq 2` and `\leq
          3`, or ``None``. This specifies the minimum degree of the generated
          quadrangulations. If this is ``None`` and the minimum connectivity is
          specified, then this is set to the same value as the minimum
          connectivity. If the minimum connectivity is also equal to ``None``,
          then this is set to 2.

        - ``minimum_connectivity`` - default: ``None`` - a value `\geq 2` and
          `\leq 3`, or ``None``. This specifies the minimum connectivity of the
          generated quadrangulations. If this is ``None`` and the option
          ``no_nonfacial_quadrangles`` is set to ``True``, then this is set to
          3. Otherwise if this is ``None`` and the minimum degree is specified,
          then this is set to the minimum degree. If the minimum degree is also
          equal to ``None``, then this is set to 3.

        - ``no_nonfacial_quadrangles`` - default: ``False`` - if ``True`` only
          quadrangulations with no non-facial quadrangles are generated. This
          option cannot be used if ``minimum_connectivity`` is set to 2.

        - ``dual`` - default: ``False`` - if ``True`` return instead the
          planar duals of the generated graphs.

        OUTPUT:

        An iterator which will produce all planar quadrangulations with the given
        number of vertices as Sage graphs with an embedding set. These will be
        simple graphs (no loops, no multiple edges, no directed edges).

        .. SEEALSO::

            - :meth:`~sage.graphs.generic_graph.GenericGraph.set_embedding`,
              :meth:`~sage.graphs.generic_graph.GenericGraph.get_embedding` --
              get/set methods for embeddings.

        EXAMPLES:

        The cube is the only 3-connected planar quadrangulation on 8 vertices::

            sage: gen = graphs.quadrangulations(8, minimum_connectivity=3)  # optional plantri
            sage: g = next(gen)                                            # optional plantri
            sage: g.is_isomorphic(graphs.CubeGraph(3))                      # optional plantri
            True
            sage: next(gen)                                                # optional plantri
            Traceback (most recent call last):
            ...
            StopIteration

        An overview of the number of quadrangulations on up to 12 vertices. This
        agrees with :oeis:`A113201`::

            sage: for i in range(4,13):                          # optional plantri
            ....:     L =  len(list(graphs.quadrangulations(i))) # optional plantri
            ....:     print("{:2d}   {:3d}".format(i,L))         # optional plantri
             4     1
             5     1
             6     2
             7     3
             8     9
             9    18
            10    62
            11   198
            12   803

        There are 2 planar quadrangulation on 12 vertices that do not have a
        non-facial quadrangle::

            sage: len([g for g in graphs.quadrangulations(12, no_nonfacial_quadrangles=True)])  # optional plantri
            2

        Setting ``dual=True`` gives the planar dual graphs::

            sage: [len(g) for g in graphs.quadrangulations(12, no_nonfacial_quadrangles=True, dual=True)]  # optional plantri
            [10, 10]
        """
        from sage.misc.package import is_package_installed
        if not is_package_installed("plantri"):
            raise TypeError("the optional plantri package is not installed")

        # number of vertices should be positive
        if order < 0:
            raise ValueError("Number of vertices should be positive.")

        # plantri can only output planar quadrangulations on up to 64 vertices
        if order > 64:
            raise ValueError("Number of vertices should be at most 64.")

        # minimum connectivity should be None, 2 or 3
        if minimum_connectivity not in {None, 2, 3}:
            raise ValueError("Minimum connectivity should be None, 2 or 3.")

        # minimum degree should be None, 2 or 3
        if minimum_degree not in {None, 2, 3}:
            raise ValueError("Minimum degree should be None, 2 or 3.")

        if (no_nonfacial_quadrangles and
            minimum_connectivity == 2):
                raise NotImplementedError("Generation of no non-facial quadrangles and minimum connectivity 2 is not implemented")

        # check combination of values of minimum degree and minimum connectivity
        if minimum_connectivity is None:
            if minimum_degree is not None:
                minimum_connectivity = min(2, minimum_degree)
            else:
                minimum_degree, minimum_connectivity = 2, 2
        else:
            if minimum_degree is None:
                minimum_degree = minimum_connectivity
            elif minimum_degree < minimum_connectivity:
                raise ValueError("Minimum connectivity can be at most the minimum degree.")

        minimum_order = {2:4, 3:8}[minimum_degree]

        if order < minimum_order:
            return

        if no_nonfacial_quadrangles:
            # for plantri -q the option -c4 means 3-connected with no non-facial quadrangles
            minimum_connectivity = 4


        cmd = 'plantri -qm{}c{}{} {}'
        command = cmd.format(minimum_degree,
                             minimum_connectivity,
                             'd' if dual else '',
                             order)

        import subprocess
        sp = subprocess.Popen(command, shell=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, close_fds=True)

        for G in graphs._read_planar_code(sp.stdout):
            yield(G)

###########################################################################
# Basic Graphs
###########################################################################
    import sage.graphs.generators.basic
    BullGraph                = staticmethod(sage.graphs.generators.basic.BullGraph)
    ButterflyGraph           = staticmethod(sage.graphs.generators.basic.ButterflyGraph)
    CircularLadderGraph      = staticmethod(sage.graphs.generators.basic.CircularLadderGraph)
    ClawGraph                = staticmethod(sage.graphs.generators.basic.ClawGraph)
    CycleGraph               = staticmethod(sage.graphs.generators.basic.CycleGraph)
    CompleteGraph            = staticmethod(sage.graphs.generators.basic.CompleteGraph)
    CompleteBipartiteGraph   = staticmethod(sage.graphs.generators.basic.CompleteBipartiteGraph)
    CompleteMultipartiteGraph= staticmethod(sage.graphs.generators.basic.CompleteMultipartiteGraph)
    DiamondGraph             = staticmethod(sage.graphs.generators.basic.DiamondGraph)
    EmptyGraph               = staticmethod(sage.graphs.generators.basic.EmptyGraph)
    Grid2dGraph              = staticmethod(sage.graphs.generators.basic.Grid2dGraph)
    GridGraph                = staticmethod(sage.graphs.generators.basic.GridGraph)
    HouseGraph               = staticmethod(sage.graphs.generators.basic.HouseGraph)
    HouseXGraph              = staticmethod(sage.graphs.generators.basic.HouseXGraph)
    LadderGraph              = staticmethod(sage.graphs.generators.basic.LadderGraph)
    LollipopGraph            = staticmethod(sage.graphs.generators.basic.LollipopGraph)
    PathGraph                = staticmethod(sage.graphs.generators.basic.PathGraph)
    StarGraph                = staticmethod(sage.graphs.generators.basic.StarGraph)
    Toroidal6RegularGrid2dGraph = staticmethod(sage.graphs.generators.basic.Toroidal6RegularGrid2dGraph)
    ToroidalGrid2dGraph      = staticmethod(sage.graphs.generators.basic.ToroidalGrid2dGraph)

###########################################################################
# Small Graphs
###########################################################################
    import sage.graphs.generators.smallgraphs
    Balaban10Cage            = staticmethod(sage.graphs.generators.smallgraphs.Balaban10Cage)
    Balaban11Cage            = staticmethod(sage.graphs.generators.smallgraphs.Balaban11Cage)
    BidiakisCube             = staticmethod(sage.graphs.generators.smallgraphs.BidiakisCube)
    BiggsSmithGraph          = staticmethod(sage.graphs.generators.smallgraphs.BiggsSmithGraph)
    BlanusaFirstSnarkGraph   = staticmethod(sage.graphs.generators.smallgraphs.BlanusaFirstSnarkGraph)
    BlanusaSecondSnarkGraph  = staticmethod(sage.graphs.generators.smallgraphs.BlanusaSecondSnarkGraph)
    BrinkmannGraph           = staticmethod(sage.graphs.generators.smallgraphs.BrinkmannGraph)
    BrouwerHaemersGraph      = staticmethod(sage.graphs.generators.smallgraphs.BrouwerHaemersGraph)
    BuckyBall                = staticmethod(sage.graphs.generators.smallgraphs.BuckyBall)
    CameronGraph             = staticmethod(sage.graphs.generators.smallgraphs.CameronGraph)
    Cell600                  = staticmethod(sage.graphs.generators.smallgraphs.Cell600)
    Cell120                  = staticmethod(sage.graphs.generators.smallgraphs.Cell120)
    ChvatalGraph             = staticmethod(sage.graphs.generators.smallgraphs.ChvatalGraph)
    ClebschGraph             = staticmethod(sage.graphs.generators.smallgraphs.ClebschGraph)
    CoxeterGraph             = staticmethod(sage.graphs.generators.smallgraphs.CoxeterGraph)
    DejterGraph              = staticmethod(sage.graphs.generators.smallgraphs.DejterGraph)
    DesarguesGraph           = staticmethod(sage.graphs.generators.smallgraphs.DesarguesGraph)
    DoubleStarSnark          = staticmethod(sage.graphs.generators.smallgraphs.DoubleStarSnark)
    DurerGraph               = staticmethod(sage.graphs.generators.smallgraphs.DurerGraph)
    DyckGraph                = staticmethod(sage.graphs.generators.smallgraphs.DyckGraph)
    EllinghamHorton54Graph   = staticmethod(sage.graphs.generators.smallgraphs.EllinghamHorton54Graph)
    EllinghamHorton78Graph   = staticmethod(sage.graphs.generators.smallgraphs.EllinghamHorton78Graph)
    ErreraGraph              = staticmethod(sage.graphs.generators.smallgraphs.ErreraGraph)
    F26AGraph                = staticmethod(sage.graphs.generators.smallgraphs.F26AGraph)
    FlowerSnark              = staticmethod(sage.graphs.generators.smallgraphs.FlowerSnark)
    FolkmanGraph             = staticmethod(sage.graphs.generators.smallgraphs.FolkmanGraph)
    FosterGraph              = staticmethod(sage.graphs.generators.smallgraphs.FosterGraph)
    FranklinGraph            = staticmethod(sage.graphs.generators.smallgraphs.FranklinGraph)
    FruchtGraph              = staticmethod(sage.graphs.generators.smallgraphs.FruchtGraph)
    GoldnerHararyGraph       = staticmethod(sage.graphs.generators.smallgraphs.GoldnerHararyGraph)
    GossetGraph              = staticmethod(sage.graphs.generators.smallgraphs.GossetGraph)
    GrayGraph                = staticmethod(sage.graphs.generators.smallgraphs.GrayGraph)
    GrotzschGraph            = staticmethod(sage.graphs.generators.smallgraphs.GrotzschGraph)
    HallJankoGraph           = staticmethod(sage.graphs.generators.smallgraphs.HallJankoGraph)
    WellsGraph               = staticmethod(sage.graphs.generators.smallgraphs.WellsGraph)
    HarborthGraph            = staticmethod(sage.graphs.generators.smallgraphs.HarborthGraph)
    HarriesGraph             = staticmethod(sage.graphs.generators.smallgraphs.HarriesGraph)
    HarriesWongGraph         = staticmethod(sage.graphs.generators.smallgraphs.HarriesWongGraph)
    HeawoodGraph             = staticmethod(sage.graphs.generators.smallgraphs.HeawoodGraph)
    HerschelGraph            = staticmethod(sage.graphs.generators.smallgraphs.HerschelGraph)
    HigmanSimsGraph          = staticmethod(sage.graphs.generators.smallgraphs.HigmanSimsGraph)
    HoffmanGraph             = staticmethod(sage.graphs.generators.smallgraphs.HoffmanGraph)
    HoffmanSingletonGraph    = staticmethod(sage.graphs.generators.smallgraphs.HoffmanSingletonGraph)
    HoltGraph                = staticmethod(sage.graphs.generators.smallgraphs.HoltGraph)
    HortonGraph              = staticmethod(sage.graphs.generators.smallgraphs.HortonGraph)
    JankoKharaghaniGraph     = staticmethod(sage.graphs.generators.smallgraphs.JankoKharaghaniGraph)
    JankoKharaghaniTonchevGraph  = staticmethod(sage.graphs.generators.smallgraphs.JankoKharaghaniTonchevGraph)
    KittellGraph             = staticmethod(sage.graphs.generators.smallgraphs.KittellGraph)
    KrackhardtKiteGraph      = staticmethod(sage.graphs.generators.smallgraphs.KrackhardtKiteGraph)
    Klein3RegularGraph       = staticmethod(sage.graphs.generators.smallgraphs.Klein3RegularGraph)
    Klein7RegularGraph       = staticmethod(sage.graphs.generators.smallgraphs.Klein7RegularGraph)
    LocalMcLaughlinGraph     = staticmethod(sage.graphs.generators.smallgraphs.LocalMcLaughlinGraph)
    LjubljanaGraph           = staticmethod(sage.graphs.generators.smallgraphs.LjubljanaGraph)
    LivingstoneGraph         = staticmethod(sage.graphs.generators.smallgraphs.LivingstoneGraph)
    M22Graph                 = staticmethod(sage.graphs.generators.smallgraphs.M22Graph)
    MarkstroemGraph          = staticmethod(sage.graphs.generators.smallgraphs.MarkstroemGraph)
    MathonStronglyRegularGraph = staticmethod(sage.graphs.generators.smallgraphs.MathonStronglyRegularGraph)
    McGeeGraph               = staticmethod(sage.graphs.generators.smallgraphs.McGeeGraph)
    McLaughlinGraph          = staticmethod(sage.graphs.generators.smallgraphs.McLaughlinGraph)
    MeredithGraph            = staticmethod(sage.graphs.generators.smallgraphs.MeredithGraph)
    MoebiusKantorGraph       = staticmethod(sage.graphs.generators.smallgraphs.MoebiusKantorGraph)
    MoserSpindle             = staticmethod(sage.graphs.generators.smallgraphs.MoserSpindle)
    NauruGraph               = staticmethod(sage.graphs.generators.smallgraphs.NauruGraph)
    PappusGraph              = staticmethod(sage.graphs.generators.smallgraphs.PappusGraph)
    PoussinGraph             = staticmethod(sage.graphs.generators.smallgraphs.PoussinGraph)
    PerkelGraph              = staticmethod(sage.graphs.generators.smallgraphs.PerkelGraph)
    PetersenGraph            = staticmethod(sage.graphs.generators.smallgraphs.PetersenGraph)
    RobertsonGraph           = staticmethod(sage.graphs.generators.smallgraphs.RobertsonGraph)
    SchlaefliGraph           = staticmethod(sage.graphs.generators.smallgraphs.SchlaefliGraph)
    ShrikhandeGraph          = staticmethod(sage.graphs.generators.smallgraphs.ShrikhandeGraph)
    SimsGewirtzGraph         = staticmethod(sage.graphs.generators.smallgraphs.SimsGewirtzGraph)
    SousselierGraph          = staticmethod(sage.graphs.generators.smallgraphs.SousselierGraph)
    SylvesterGraph           = staticmethod(sage.graphs.generators.smallgraphs.SylvesterGraph)
    SzekeresSnarkGraph       = staticmethod(sage.graphs.generators.smallgraphs.SzekeresSnarkGraph)
    ThomsenGraph             = staticmethod(sage.graphs.generators.smallgraphs.ThomsenGraph)
    TietzeGraph              = staticmethod(sage.graphs.generators.smallgraphs.TietzeGraph)
    Tutte12Cage              = staticmethod(sage.graphs.generators.smallgraphs.Tutte12Cage)
    TruncatedIcosidodecahedralGraph = staticmethod(sage.graphs.generators.smallgraphs.TruncatedIcosidodecahedralGraph)
    TruncatedTetrahedralGraph= staticmethod(sage.graphs.generators.smallgraphs.TruncatedTetrahedralGraph)
    TutteCoxeterGraph        = staticmethod(sage.graphs.generators.smallgraphs.TutteCoxeterGraph)
    TutteGraph               = staticmethod(sage.graphs.generators.smallgraphs.TutteGraph)
    WagnerGraph              = staticmethod(sage.graphs.generators.smallgraphs.WagnerGraph)
    WatkinsSnarkGraph        = staticmethod(sage.graphs.generators.smallgraphs.WatkinsSnarkGraph)
    WienerArayaGraph         = staticmethod(sage.graphs.generators.smallgraphs.WienerArayaGraph)
    SuzukiGraph              = staticmethod(sage.graphs.generators.smallgraphs.SuzukiGraph)

###########################################################################
# Platonic Solids
###########################################################################
    import sage.graphs.generators.platonic_solids
    DodecahedralGraph        = staticmethod(sage.graphs.generators.platonic_solids.DodecahedralGraph)
    HexahedralGraph          = staticmethod(sage.graphs.generators.platonic_solids.HexahedralGraph)
    IcosahedralGraph         = staticmethod(sage.graphs.generators.platonic_solids.IcosahedralGraph)
    OctahedralGraph          = staticmethod(sage.graphs.generators.platonic_solids.OctahedralGraph)
    TetrahedralGraph         = staticmethod(sage.graphs.generators.platonic_solids.TetrahedralGraph)

###########################################################################
# Families
###########################################################################
    import sage.graphs.generators.families
    BalancedTree           = staticmethod(sage.graphs.generators.families.BalancedTree)
    BarbellGraph           = staticmethod(sage.graphs.generators.families.BarbellGraph)
    BubbleSortGraph        = staticmethod(sage.graphs.generators.families.BubbleSortGraph)
    chang_graphs           = staticmethod(sage.graphs.generators.families.chang_graphs)
    CirculantGraph         = staticmethod(sage.graphs.generators.families.CirculantGraph)
    CubeGraph              = staticmethod(sage.graphs.generators.families.CubeGraph)
    DorogovtsevGoltsevMendesGraph = staticmethod(sage.graphs.generators.families.DorogovtsevGoltsevMendesGraph)
    FibonacciTree          = staticmethod(sage.graphs.generators.families.FibonacciTree)
    FoldedCubeGraph        = staticmethod(sage.graphs.generators.families.FoldedCubeGraph)
    FriendshipGraph        = staticmethod(sage.graphs.generators.families.FriendshipGraph)
    FuzzyBallGraph         = staticmethod(sage.graphs.generators.families.FuzzyBallGraph)
    GeneralizedPetersenGraph = staticmethod(sage.graphs.generators.families.GeneralizedPetersenGraph)
    GoethalsSeidelGraph    = staticmethod(sage.graphs.generators.families.GoethalsSeidelGraph)
    HanoiTowerGraph        = staticmethod(sage.graphs.generators.families.HanoiTowerGraph)
    HararyGraph            = staticmethod(sage.graphs.generators.families.HararyGraph)
    HyperStarGraph         = staticmethod(sage.graphs.generators.families.HyperStarGraph)
    JohnsonGraph           = staticmethod(sage.graphs.generators.families.JohnsonGraph)
    KneserGraph            = staticmethod(sage.graphs.generators.families.KneserGraph)
    LCFGraph               = staticmethod(sage.graphs.generators.families.LCFGraph)
    line_graph_forbidden_subgraphs = staticmethod(sage.graphs.generators.families.line_graph_forbidden_subgraphs)
    MathonPseudocyclicMergingGraph = staticmethod(sage.graphs.generators.families.MathonPseudocyclicMergingGraph)
    MathonPseudocyclicStronglyRegularGraph = staticmethod(sage.graphs.generators.families.MathonPseudocyclicStronglyRegularGraph)
    MycielskiGraph         = staticmethod(sage.graphs.generators.families.MycielskiGraph)
    MycielskiStep          = staticmethod(sage.graphs.generators.families.MycielskiStep)
    NKStarGraph            = staticmethod(sage.graphs.generators.families.NKStarGraph)
    NStarGraph             = staticmethod(sage.graphs.generators.families.NStarGraph)
    OddGraph               = staticmethod(sage.graphs.generators.families.OddGraph)
    PaleyGraph             = staticmethod(sage.graphs.generators.families.PaleyGraph)
    PasechnikGraph         = staticmethod(sage.graphs.generators.families.PasechnikGraph)
    petersen_family        = staticmethod(sage.graphs.generators.families.petersen_family)
    RingedTree             = staticmethod(sage.graphs.generators.families.RingedTree)
    SierpinskiGasketGraph  = staticmethod(sage.graphs.generators.families.SierpinskiGasketGraph)
    SquaredSkewHadamardMatrixGraph = staticmethod(sage.graphs.generators.families.SquaredSkewHadamardMatrixGraph)
    SwitchedSquaredSkewHadamardMatrixGraph = staticmethod(sage.graphs.generators.families.SwitchedSquaredSkewHadamardMatrixGraph)
    strongly_regular_graph = staticmethod(sage.graphs.strongly_regular_db.strongly_regular_graph)
    trees                  = staticmethod(sage.graphs.generators.families.trees)
    TuranGraph             = staticmethod(sage.graphs.generators.families.TuranGraph)
    WheelGraph             = staticmethod(sage.graphs.generators.families.WheelGraph)

###########################################################################
# Graphs from classical geometries over `F_q`
###########################################################################
    import sage.graphs.generators.classical_geometries
    AffineOrthogonalPolarGraph = staticmethod(sage.graphs.generators.classical_geometries.AffineOrthogonalPolarGraph)
    AhrensSzekeresGeneralizedQuadrangleGraph = staticmethod(sage.graphs.generators.classical_geometries.AhrensSzekeresGeneralizedQuadrangleGraph)
    NonisotropicOrthogonalPolarGraph = staticmethod(sage.graphs.generators.classical_geometries.NonisotropicOrthogonalPolarGraph)
    NonisotropicUnitaryPolarGraph = staticmethod(sage.graphs.generators.classical_geometries.NonisotropicUnitaryPolarGraph)
    OrthogonalPolarGraph   = staticmethod(sage.graphs.generators.classical_geometries.OrthogonalPolarGraph)
    SymplecticDualPolarGraph = staticmethod(sage.graphs.generators.classical_geometries.SymplecticDualPolarGraph)
    SymplecticGraph   = staticmethod(sage.graphs.generators.classical_geometries.SymplecticGraph)
    SymplecticPolarGraph   = staticmethod(sage.graphs.generators.classical_geometries.SymplecticPolarGraph)
    TaylorTwographDescendantSRG = \
             staticmethod(sage.graphs.generators.classical_geometries.TaylorTwographDescendantSRG)
    TaylorTwographSRG      = staticmethod(sage.graphs.generators.classical_geometries.TaylorTwographSRG)
    T2starGeneralizedQuadrangleGraph      = staticmethod(sage.graphs.generators.classical_geometries.T2starGeneralizedQuadrangleGraph)
    HaemersGraph      = staticmethod(sage.graphs.generators.classical_geometries.HaemersGraph)
    CossidentePenttilaGraph = staticmethod(sage.graphs.generators.classical_geometries.CossidentePenttilaGraph)
    UnitaryDualPolarGraph  = staticmethod(sage.graphs.generators.classical_geometries.UnitaryDualPolarGraph)
    UnitaryPolarGraph      = staticmethod(sage.graphs.generators.classical_geometries.UnitaryPolarGraph)

###########################################################################
# Chessboard Graphs
###########################################################################
    import sage.graphs.generators.chessboard
    ChessboardGraphGenerator = staticmethod(sage.graphs.generators.chessboard.ChessboardGraphGenerator)
    BishopGraph              = staticmethod(sage.graphs.generators.chessboard.BishopGraph)
    KingGraph                = staticmethod(sage.graphs.generators.chessboard.KingGraph)
    KnightGraph              = staticmethod(sage.graphs.generators.chessboard.KnightGraph)
    QueenGraph               = staticmethod(sage.graphs.generators.chessboard.QueenGraph)
    RookGraph                = staticmethod(sage.graphs.generators.chessboard.RookGraph)

###########################################################################
# Intersection graphs
###########################################################################
    import sage.graphs.generators.intersection
    IntervalGraph            = staticmethod(sage.graphs.generators.intersection.IntervalGraph)
    IntersectionGraph        = staticmethod(sage.graphs.generators.intersection.IntersectionGraph)
    PermutationGraph         = staticmethod(sage.graphs.generators.intersection.PermutationGraph)
    OrthogonalArrayBlockGraph  = staticmethod(sage.graphs.generators.intersection.OrthogonalArrayBlockGraph)
    ToleranceGraph           = staticmethod(sage.graphs.generators.intersection.ToleranceGraph)

###########################################################################
# Random Graphs
###########################################################################
    import sage.graphs.generators.random
    RandomBarabasiAlbert     = staticmethod(sage.graphs.generators.random.RandomBarabasiAlbert)
    RandomBipartite          = staticmethod(sage.graphs.generators.random.RandomBipartite)
    RandomBoundedToleranceGraph = staticmethod(sage.graphs.generators.random.RandomBoundedToleranceGraph)
    RandomGNM                = staticmethod(sage.graphs.generators.random.RandomGNM)
    RandomGNP                = staticmethod(sage.graphs.generators.random.RandomGNP)
    RandomHolmeKim           = staticmethod(sage.graphs.generators.random.RandomHolmeKim)
    RandomIntervalGraph      = staticmethod(sage.graphs.generators.random.RandomIntervalGraph)
    RandomLobster            = staticmethod(sage.graphs.generators.random.RandomLobster)
    RandomNewmanWattsStrogatz = staticmethod(sage.graphs.generators.random.RandomNewmanWattsStrogatz)
    RandomRegular            = staticmethod(sage.graphs.generators.random.RandomRegular)
    RandomShell              = staticmethod(sage.graphs.generators.random.RandomShell)
    RandomToleranceGraph     = staticmethod(sage.graphs.generators.random.RandomToleranceGraph)
    RandomTreePowerlaw       = staticmethod(sage.graphs.generators.random.RandomTreePowerlaw)
    RandomTree               = staticmethod(sage.graphs.generators.random.RandomTree)
    RandomTriangulation      = staticmethod(sage.graphs.generators.random.RandomTriangulation)

###########################################################################
# World Map
###########################################################################
    import sage.graphs.generators.world_map
    WorldMap = staticmethod(sage.graphs.generators.world_map.WorldMap)

###########################################################################
# Degree Sequence
###########################################################################
    import sage.graphs.generators.degree_sequence
    DegreeSequence           = staticmethod(sage.graphs.generators.degree_sequence.DegreeSequence)
    DegreeSequenceBipartite  = staticmethod(sage.graphs.generators.degree_sequence.DegreeSequenceBipartite)
    DegreeSequenceConfigurationModel = staticmethod(sage.graphs.generators.degree_sequence.DegreeSequenceConfigurationModel)
    DegreeSequenceTree       = staticmethod(sage.graphs.generators.degree_sequence.DegreeSequenceTree)
    DegreeSequenceExpected   = staticmethod(sage.graphs.generators.degree_sequence.DegreeSequenceExpected)

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
                        if m_z.relabel(dict(enumerate(possibility)), check_input=False, inplace=False) == g:
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
