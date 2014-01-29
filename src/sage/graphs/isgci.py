r"""
Information System on Graph Classes and their Inclusions

This module implements an interface to the `ISGCI <http://www.graphclasses.org/>`_ database in Sage.

This database gathers information on graph classes and their
inclusions in each other. It also contains information on the
complexity of several computational problems.

It is available on the `GraphClasses.org <http://www.graphclasses.org/>`_
website maintained by H.N. de Ridder et al.

How to use it?
--------------

The current limited aim of this module is to provide a pure Sage/Python
interface which is easier to deal with than a XML file. I hope it will be
rewritten many times until it stabilizes :-)

Presently, it is possible to use this database through the variables and methods
present in the :obj:`graph_classes <GraphClasses>` object.
For instance::

    sage: Trees = graph_classes.Tree
    sage: Chordal = graph_classes.Chordal

It is then possible to check the inclusion of classes inside of others, if the
information is available in the database::

    sage: Trees <= Chordal
    True

And indeed, trees are chordal graphs.

.. WARNING::

    The ISGCI database is not all-knowing, and so comparing two classes can
    return ``True``, ``False``, or ``Unknown`` (see the :mod:`documentation of
    the Unknown truth value <sage.misc.unknown>`).

    An *unknown* answer to ``A <= B`` only means that ISGCI cannot deduce from
    the information in its database that ``A`` is a subclass of ``B`` nor that
    it is not. For instance, ISGCI does not know at the moment that some chordal
    graphs are not trees::

        sage: graph_classes.Chordal <= graph_classes.Tree
        Unknown

Given a graph class, one can obtain its associated information in the
ISGCI database with the :meth:`~GraphClass.description` method::

    sage: Chordal.description()
    Class of graphs : Chordal
    -------------------------
    type                           :  base
    ID                             :  gc_32
    name                           :  chordal
    <BLANKLINE>
    Problems :
    -----------
    3-Colourability                :  Linear
    Clique                         :  Polynomial
    Clique cover                   :  Polynomial
    Cliquewidth                    :  Unbounded
    Cliquewidth expression         :  NP-complete
    Colourability                  :  Linear
    Domination                     :  NP-complete
    Independent set                :  Linear
    Recognition                    :  Linear
    Treewidth                      :  Polynomial
    Weighted clique                :  Polynomial
    Weighted independent set       :  Linear

It is possible to obtain the complete list of the classes stored in ISGCI by
calling the :meth:`~GraphClasses.show_all` method (beware -- long output)::

    sage: graph_classes.show_all()
    ID        | name                                     | type                 | smallgraph
    ----------------------------------------------------------------------------------------------------------------------
    gc_309    | $K_4$--minor--free                       | base                 |
    gc_541    | $N^*$                                    | base                 |
    gc_215    | $N^*$--perfect                           | base                 |
    gc_5      | $P_4$--bipartite                         | base                 |
    gc_3      | $P_4$--brittle                           | base                 |
    gc_6      | $P_4$--comparability                     | base                 |
    gc_7      | $P_4$--extendible                        | base                 |
    ...

Until a proper search method is implemented, this lets one find
classes which do not appear in :obj:`graph_classes.* <GraphClasses>`.

To retrieve a class of graph from its ISGCI ID one may use
the :meth:`~GraphClasses.get_class` method::

    sage: GC = graph_classes.get_class("gc_5")
    sage: GC
    $P_4$--bipartite graphs

Predefined classes
------------------

:obj:`graph_classes <GraphClasses>` currently predefines the following graph classes

.. list-table::
   :widths: 20 30
   :header-rows: 1

   * - Class
     - Related methods

   * - BinaryTrees

     - :meth:`~sage.graphs.graph_generators.GraphGenerators.BalancedTree`,
       :meth:`~Graph.is_tree`

   * - Bipartite

     - :meth:`~sage.graphs.graph_generators.GraphGenerators.BalancedTree`,
       :meth:`~sage.graphs.graph.Graph.is_bipartite`

   * - Block

     - :meth:`~sage.graphs.generic_graph.GenericGraph.blocks_and_cut_vertices`,

   * - Chordal

     - :meth:`~sage.graphs.generic_graph.GenericGraph.is_chordal`

   * - Comparability
     -

   * - Gallai

     - :meth:`~sage.graphs.generic_graph.GenericGraph.is_gallai_tree`

   * - Grid

     - :meth:`~sage.graphs.graph_generators.GraphGenerators.Grid2dGraph`,
       :meth:`~sage.graphs.graph_generators.GraphGenerators.GridGraph`

   * - Interval

     - :meth:`~sage.graphs.graph_generators.GraphGenerators.RandomInterval`,
       :meth:`~sage.graphs.graph_generators.GraphGenerators.IntervalGraph`,
       :meth:`~sage.graphs.generic_graph.GenericGraph.is_interval`

   * - Line

     - :meth:`~sage.graphs.graph_generators.GraphGenerators.line_graph_forbidden_subgraphs`,
       :meth:`~sage.graphs.graph.Graph.is_line_graph`

   * - Modular

     - :meth:`~sage.graphs.graph.Graph.modular_decomposition`

   * - Outerplanar

     - :meth:`~sage.graphs.generic_graph.GenericGraph.is_circular_planar`

   * - Perfect

     - :meth:`~sage.graphs.graph.Graph.is_perfect`

   * - Planar

     - :meth:`~sage.graphs.generic_graph.GenericGraph.is_planar`

   * - Split

     - :meth:`~sage.graphs.graph.Graph.is_split`

   * - Tree

     - :meth:`~sage.graphs.graph_generators.GraphGenerators.trees`,
       :meth:`~Graph.is_tree`

   * - UnitDisk
     - :meth:`~sage.graphs.graph_generators.GraphGenerators.IntervalGraph`,

   * - UnitInterval

     - :meth:`~sage.graphs.generic_graph.GenericGraph.is_interval`

Sage's view of ISGCI
--------------------

The database is stored by Sage in two ways.

**The classes**: the list of all graph classes and their properties is stored
in a huge dictionary (see :meth:`~sage.graphs.isgci.GraphClasses.classes`).
Below is what Sage knows of ``gc_249``::

    sage: graph_classes.classes()['gc_249']        # random
    {'problems':
        {'Independent set': 'Polynomial',
         'Treewidth': 'Unknown',
         'Weighted independent set': 'Polynomial',
         'Cliquewidth expression': 'NP-complete',
         'Weighted clique': 'Polynomial',
         'Clique cover': 'Unknown',
         'Domination': 'NP-complete',
         'Clique': 'Polynomial',
         'Colourability': 'NP-complete',
         'Cliquewidth': 'Unbounded',
         '3-Colourability': 'NP-complete',
         'Recognition': 'Linear'},
     'type': 'base',
     'ID': 'gc_249',
     'name': 'line'}

**The class inclusion digraph**: Sage remembers the class inclusions through
the inclusion digraph (see :meth:`~sage.graphs.isgci.GraphClasses.inclusion_digraph`).
Its nodes are ID of ISGCI classes::

    sage: d = graph_classes.inclusion_digraph()
    sage: d.vertices()[-10:]
    ['gc_990', 'gc_991', 'gc_992', 'gc_993', 'gc_994', 'gc_995', 'gc_996', 'gc_997', 'gc_998', 'gc_999']

An arc from ``gc1`` to ``gc2`` means that ``gc1`` is a superclass of
``gc2``. This being said, not all edges are stored ! To ensure that a given
class is included in another one, we have to check whether there is in the
digraph a ``path`` from the first one to the other::

    sage: bip_id = graph_classes.Bipartite._gc_id
    sage: perfect_id = graph_classes.Perfect._gc_id
    sage: d.has_edge(perfect_id, bip_id)
    False
    sage: d.distance(perfect_id, bip_id)
    2

Hence bipartite graphs are perfect graphs. We can see how ISGCI obtains this
result ::

    sage: p = d.shortest_path(perfect_id, bip_id)
    sage: len(p) - 1
    2
    sage: print p                  # random
    ['gc_56', 'gc_76', 'gc_69']
    sage: for c in p:
    ...      print graph_classes.get_class(c)
    perfect graphs
    ...
    bipartite graphs

What ISGCI knows is that perfect graphs contain unimodular graph which contain
bipartite graphs. Therefore bipartite graphs are perfect !

.. note::

    The inclusion digraph is **NOT ACYCLIC**. Indeed, several entries
    exist in the ISGCI database which represent the same graph class,
    for instance Perfect graphs and Berge graphs::

        sage: graph_classes.inclusion_digraph().is_directed_acyclic()
        False
        sage: Berge = graph_classes.get_class("gc_274"); Berge
        Berge graphs
        sage: Perfect = graph_classes.get_class("gc_56"); Perfect
        perfect graphs
        sage: Berge <= Perfect
        True
        sage: Perfect <= Berge
        True
        sage: Perfect == Berge
        True

Information for developpers
----------------------------

* The database is loaded not *so* large, but it is still preferable to
  only load it on demand. This is achieved through the cached methods
  :meth:`~sage.graphs.isgci.GraphClasses.classes` and
  :meth:`~sage.graphs.isgci.GraphClasses.inclusion_digraph`.

* Upon the first access to the database, the information is extracted
  from the XML file and stored in the cache of three methods:

  * ``sage.graphs.isgci._classes`` (dictionary)
  * ``sage.graphs.isgci._inclusions`` (list of dictionaries)
  * ``sage.graphs.isgci._inclusion_digraph`` (DiGraph)

  Note that the digraph is only built if necessary (for instance if
  the user tries to compare two classes).

.. todo::

    Technical things:

    * Query the database for non-inclusion results so that comparisons can
      return ``False``, and implement strict inclusions.

    * Implement a proper search method for the classes not listed in
      :obj:`graph_classes <GraphClasses>`

      .. seealso: :func:`sage.graphs.isgci.show_all`.

    * Some of the graph classes appearing in :obj:`graph_classes
      <GraphClasses>` already have a recognition
      algorithm implemented in Sage. It would be so nice to be able to
      write ``g in Trees``, ``g in Perfect``, ``g in Chordal``, ... :-)

    Long-term stuff:

    * Implement simple accessors for all the information in the ISGCI
      database (as can be done from the website)

    * Implement intersection of graph classes

    * Write generic recognition algorithms for specific classes (when a graph class
      is defined by the exclusion of subgraphs, one can write a generic algorithm
      checking the existence of each of the graphs, and this method already exists
      in Sage).

    * Improve the performance of Sage's graph library by letting it
      take advantage of the properties of graph classes. For example,
      :meth:`Graph.independent_set` could use the library to detect
      that a given graph is, say, a tree or a planar graph, and use a
      specialized algorithm for finding an independent set.

AUTHORS:
--------

* H.N. de Ridder et al. (ISGCI database)
* Nathann Cohen (Sage implementation)

Methods
-------
"""

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import CachedRepresentation, UniqueRepresentation
from sage.misc.unknown import Unknown
from sage.misc.misc import SAGE_SHARE

#*****************************************************************************
#      Copyright (C) 2011 Nathann Cohen <nathann.cohen@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

_XML_FILE = "isgci_sage.xml"

class GraphClass(SageObject, CachedRepresentation):
    r"""
    An instance of this class represents a Graph Class, matching some entry in
    the ISGCI database.

    EXAMPLE:

    Testing the inclusion of two classes::

        sage: Chordal = graph_classes.Chordal
        sage: Trees = graph_classes.Tree
        sage: Trees <= Chordal
        True
        sage: Chordal <= Trees
        Unknown

    TEST::

        sage: Trees >= Chordal
        Unknown
        sage: Chordal >= Trees
        True
    """
    def __init__(self, name, gc_id):
        r"""
        Class constructor

        EXAMPLE::

            sage: graph_classes.Chordal  # indirect doctest
            Chordal graphs
        """
        self._name = name
        self._gc_id = gc_id

    def _repr_(self):
        r"""
        Returns a short description of the class

        EXAMPLE::

            sage: graph_classes.Chordal  # indirect doctest
            Chordal graphs
        """
        return self._name+" graphs"

    def __hash__(self):
        r"""
        Returns the class' ID hash

        EXAMPLE::

            sage: hash(graph_classes.Chordal) == hash(graph_classes.Chordal)
            True
        """
        return hash(self._gc_id)

    def __le__(self, other):
        r"""
        <= operator

        EXAMPLE::

            sage: graph_classes.Chordal <= graph_classes.Tree
            Unknown
        """
        return other.__ge__(self)

    def __ge__(self, other):
        r"""
        >= operator

        EXAMPLE::

            sage: graph_classes.Chordal >= graph_classes.Tree
            True
        """

        inclusion_digraph = GraphClasses().inclusion_digraph()
        if inclusion_digraph.shortest_path(self._gc_id,other._gc_id) != []:
            return True
        else:
            return Unknown

    def __eq__(self, other):
        r"""
        == operator

        EXAMPLE::

            sage: graph_classes.Chordal == graph_classes.Tree
            Unknown
        """
        return self.__ge__(other) and other.__ge__(self)

    def __lt__(self, other):
        r"""
        >, !=, and < operators

        EXAMPLE::

            sage: graph_classes.Chordal > graph_classes.Tree
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: graph_classes.Chordal < graph_classes.Tree
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: graph_classes.Chordal != graph_classes.Tree
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    __gt__ = __ne__ = __lt__

    def description(self):
        r"""
        Prints the information of ISGCI about the current class.

        EXAMPLE::

            sage: graph_classes.Chordal.description()
            Class of graphs : Chordal
            -------------------------
            type                           :  base
            ID                             :  gc_32
            name                           :  chordal
            <BLANKLINE>
            Problems :
            -----------
            3-Colourability                :  Linear
            Clique                         :  Polynomial
            Clique cover                   :  Polynomial
            Cliquewidth                    :  Unbounded
            Cliquewidth expression         :  NP-complete
            Colourability                  :  Linear
            Domination                     :  NP-complete
            Independent set                :  Linear
            Recognition                    :  Linear
            Treewidth                      :  Polynomial
            Weighted clique                :  Polynomial
            Weighted independent set       :  Linear
        """

        classes = GraphClasses().classes()

        cls = classes[self._gc_id]
        print "Class of graphs : "+self._name
        print "-"*(len(self._name)+18)

        for key, value in cls.iteritems():
            if value != "" and key != "problems":
                print "{0:30} : ".format(key),
                print value

        print "\nProblems :"
        print "-"*11
        for key, value in sorted(cls["problems"].iteritems()):
            if value != "":
                print "{0:30} : ".format(key),
                print value

from sage.misc.cachefunc import cached_method

class GraphClasses(UniqueRepresentation):
    def get_class(self, id):
        r"""
        Returns the class corresponding to the given id in the ISGCI database.

        INPUT:

        - ``id`` (string) -- the desired class' ID

        .. seealso:

            :meth:`~sage.graphs.isgci.GraphClasses.show_all`

        EXAMPLE:

        With an existing id::

            sage: Cographs = graph_classes.get_class("gc_151")
            sage: Cographs
            cograph graphs

        With a wrong id::

            sage: graph_classes.get_class(-1)
            Traceback (most recent call last):
            ...
            ValueError: The given class id does not exist in the ISGCI database. Is the db too old ? You can update it with graph_classes.update_db().
        """
        classes = self.classes()
        if id in classes:
            c = classes[id]

            if "name" in c and c["name"] != "":
                name = c["name"]
            else:
                name = "class "+str(id)

            return GraphClass(name, id)
        else:
            raise ValueError("The given class id does not exist in the ISGCI database. Is the db too old ? You can update it with graph_classes.update_db().")

    @cached_method
    def classes(self):
        r"""
        Returns the graph classes, as a dictionary.

        Upon the first call, this loads the database from the local
        XML file. Subsequent calls are cached.

        EXAMPLES::

            sage: t = graph_classes.classes()
            sage: type(t)
            <type 'dict'>
            sage: sorted(t["gc_151"].keys())
            ['ID', 'name', 'problems', 'type']
            sage: t["gc_151"]['name']
            'cograph'
            sage: t["gc_151"]['problems']['Clique']
            'Linear'
        """
        classes, inclusions = self._get_ISGCI()
        self.inclusions.set_cache(inclusions)
        return classes

    @cached_method
    def inclusions(self):
        r"""
        Returns the graph class inclusions

        OUTPUT:

        a list of dictionaries

        Upon the first call, this loads the database from the local
        XML file. Subsequent calls are cached.

        EXAMPLES::

            sage: t = graph_classes.inclusions()
            sage: type(t)
            <type 'list'>
            sage: t[0]
            {'super': 'gc_2', 'sub': 'gc_1'}
        """
        self.classes()
        return self.inclusions()

    @cached_method
    def inclusion_digraph(self):
        r"""
        Returns the class inclusion digraph

        Upon the first call, this loads the database from the local
        XML file. Subsequent calls are cached.

        EXAMPLES::

            sage: g = graph_classes.inclusion_digraph(); g
            Digraph on ... vertices
        """
        classes    = self.classes()
        inclusions = self.inclusions()

        from sage.graphs.digraph import DiGraph
        inclusion_digraph = DiGraph()
        inclusion_digraph.add_vertices(classes.keys())

        for edge in inclusions:
            if edge.get("confidence","") == "unpublished":
                continue
            inclusion_digraph.add_edge(edge['super'], edge['sub'])

        return inclusion_digraph

    def _download_db(self):
        r"""
        Downloads the current version of the ISGCI db

        EXAMPLE::

            sage: graph_classes._download_db() # Not tested -- requires internet
        """

        from sage.misc.misc import SAGE_TMP
        import urllib2
        import os.path
        u = urllib2.urlopen('http://www.graphclasses.org/data.zip')
        localFile = open(os.path.join(SAGE_TMP,'isgci.zip'), 'w')
        localFile.write(u.read())
        localFile.close()
        import os, zipfile
        z = zipfile.ZipFile(os.path.join(SAGE_TMP,'isgci.zip'))

        # Save a systemwide updated copy whenever possible

        try:
            z.extract(_XML_FILE, os.path.join(SAGE_SHARE,'graphs'))
        except IOError:
            z.extract(_XML_FILE, SAGE_TMP)


    def _parse_db(self, xml_file):
        r"""
        Parses the ISGCI database and returns its content as Python objects.

        INPUT:

        - ``xml_file`` -- the name of an XML file containing the data from ISGCI

        EXAMPLE::

            sage: import os
            sage: from sage.misc.misc import SAGE_SHARE
            sage: map(type, graph_classes._parse_db(os.path.join(SAGE_SHARE,'graphs','isgci_sage.xml')))
            [<type 'dict'>, <type 'list'>]
        """
        import xml.dom.minidom
        from xml.dom.minidom import Node

        # This method is there to parse the XML file containing the ISGCI
        # database. It is admittedly not very pretty, but it builds the class we
        # want from the XML file and that's more or less all we ask it to do :-p

        doc = xml.dom.minidom.parse(xml_file)

        classes = {}

        giveme = lambda x,y : str(x.getAttribute(y))

        for node in doc.getElementsByTagName("GraphClass"):
            ID = str(node.getAttribute("id"))
            Dl = {}
            smallgraph = []
            Dl["ID"] = giveme(node, "id")
            problems = {}
            for node2 in node.childNodes:
                name = str(node2.nodeName)
                if name == "name":
                    Dl[name] = str(node2.childNodes[0].nodeValue)

                elif name ==  "smallgraph":
                    smallgraph.append(str(node2.childNodes[0].nodeValue))

                elif name == "problem":
                    problems[giveme(node2, "name")] = giveme(node2, "complexity")

            Dl["problems"] = problems

            if smallgraph:
                Dl["smallgraph"] = smallgraph

            if giveme(node, "type"):
                Dl["type"] = giveme(node, "type")


            classes[giveme(node, "id")] = Dl

        inclusions = []
        for node in doc.getElementsByTagName("incl"):
            Dl = {}
            for name in ["proper", "confidence", "super", "sub"]:
                if giveme(node, name):
                    Dl[name] = giveme(node, name)

            for node2 in node.childNodes:
                name = str(node2.nodeName)
                if name == "ref":
                    Dl[name] = str(node2.childNodes[0].nodeValue)

            inclusions.append(Dl)

        return classes, inclusions

    def update_db(self):
        r"""
        Updates the ISGCI database by downloading the latest version from internet.

        This method downloads the ISGCI database from the website
        `GraphClasses.org <http://www.graphclasses.org/>`_. It then extracts the
        zip file and parses its XML content.

        Depending on the credentials of the user running Sage when this command
        is run, one attempt is made at saving the result in Sage's directory so
        that all users can benefit from it. If the credentials are not
        sufficient, the XML file are saved instead in the user's directory (in
        the SAGE_DB folder).

        EXAMPLE::

            sage: graph_classes.update_db() # Not tested -- requires internet
        """
        from sage.misc.misc import SAGE_TMP, SAGE_DB

        self._download_db()

        print "Database downloaded"

        self.classes.clear_cache()
        self.inclusions.clear_cache()
        self.inclusion_digraph.clear_cache()

    def _get_ISGCI(self):
        r"""
        Returns the contents of the ISGCI database.

        This method is mostly for internal use, but often provides useful
        information during debugging operations.

        OUTPUT:

        A pair ``(classes, inclusions)`` where ``classes`` is a dict of dict, and
        ``inclusions`` is a list of dicts.

        .. NOTE::

            This method returns the data contained in the most recent ISGCI database
            present on the computer. See :meth:`update_db` to update the latter.

        EXAMPLE::

            sage: classes, inclusions = graph_classes._get_ISGCI()  # long time (4s on sage.math, 2012)
        """

        import os.path
        from sage.all import save, load
        from sage.misc.misc import SAGE_TMP, SAGE_DB

        try:
            open(os.path.join(SAGE_DB,_XML_FILE))

            # Which copy is the most recent on the disk ?
            if (os.path.getmtime(os.path.join(SAGE_DB,_XML_FILE)) >
                os.path.getmtime(os.path.join(SAGE_SHARE,'graphs',_XML_FILE))):

                filename = os.path.join(SAGE_DB,_XML_FILE)

            else:
                filename = os.path.join(SAGE_SHARE,'graphs',_XML_FILE)

        except IOError as e:
            filename = os.path.join(SAGE_SHARE,'graphs',_XML_FILE)

        return self._parse_db(filename)

    def show_all(self):
        r"""
        Prints all graph classes stored in ISGCI

        EXAMPLE::

            sage: graph_classes.show_all()
            ID        | name                                     | type                 | smallgraph
            ----------------------------------------------------------------------------------------------------------------------
            gc_309    | $K_4$--minor--free                       | base                 |
            gc_541    | $N^*$                                    | base                 |
            gc_215    | $N^*$--perfect                           | base                 |
            gc_5      | $P_4$--bipartite                         | base                 |
            gc_3      | $P_4$--brittle                           | base                 |
            gc_6      | $P_4$--comparability                     | base                 |
            gc_7      | $P_4$--extendible                        | base                 |
            ...
        """
        classes = self.classes()
        classes_list = classes.values()

        # We want to print the different fields, and this dictionary stores the
        # maximal number of characters of each field.
        MAX = {
            "ID" : 0,
            "type" : 0,
            "smallgraph": 0,
            "name": 0
            }

        # We sort the classes alphabetically, though we would like to display the
        # meaningful classes at the top of the list
        classes_list.sort(key = lambda x:x.get("name","zzzzz")+"{0:4}".format(int(x["ID"].split('_')[1])))

        # Maximum width of a field
        MAX_LEN = 40

        # Computing te max of each field with the database
        for key in MAX:
            MAX[key] = len(max(map(lambda x:str(x.get(key,"")),classes_list), key = len))

        # At most MAX characters per field
        for key, length in MAX.iteritems():
            MAX[key] = min(length, MAX_LEN)

        # Head of the table
        print ("{0:"+str(MAX["ID"])+"} | {1:"+str(MAX["name"])+"} | {2:"+str(MAX["type"])+"} | {3:"+str(MAX["smallgraph"])+"}").format("ID", "name", "type", "smallgraph")
        print "-"*(sum(MAX.values())+9)

        # Entries
        for entry in classes_list:
            ID = entry.get("ID","")
            name = entry.get("name","")
            type = entry.get("type","")
            smallgraph = entry.get("smallgraph","")
            print ("{0:"+str(MAX["ID"])+"} | {1:"+str(MAX["name"])+"} | {2:"+str(MAX["type"])+"} | ").format(ID, name[:MAX_LEN], type[:MAX_LEN])+str(smallgraph)[:MAX_LEN]

graph_classes = GraphClasses()

# Any object added to this list should also appear in the class' documentation, at the top of the file.
graph_classes.BinaryTrees = GraphClass("BinaryTrees", "gc_847")
graph_classes.Bipartite = GraphClass("Bipartite", "gc_69")
graph_classes.Block = GraphClass("Block", "gc_93")
graph_classes.Chordal = GraphClass("Chordal", "gc_32")
graph_classes.Comparability = GraphClass("Comparability", "gc_72")
graph_classes.Gallai = GraphClass("Gallai", "gc_73")
graph_classes.Grid = GraphClass("Grid", "gc_464")
graph_classes.Interval = GraphClass("Interval", "gc_234")
graph_classes.Line = GraphClass("Line", "gc_249")
graph_classes.Modular = GraphClass("Modular", "gc_50")
graph_classes.Outerplanar = GraphClass("Outerplanar", "gc_110")
graph_classes.Perfect = GraphClass("Perfect", "gc_56")
graph_classes.Planar = GraphClass("Planar", "gc_43")
graph_classes.Split = GraphClass("Split", "gc_39")
graph_classes.Tree = GraphClass("Tree", "gc_342")
graph_classes.UnitDisk = GraphClass("UnitDisk", "gc_389")
graph_classes.UnitInterval = GraphClass("UnitInterval", "gc_299")
