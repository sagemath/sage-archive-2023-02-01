r"""
Graph Database Module

INFO:

    This module implements classes (GraphDatabase, GraphQuery, GenericGraphQuery)
    for interfacing with the sqlite database graphs.db.

    The GraphDatabase class interfaces with the sqlite database graphs.db.
    It is an immutable database that inherits from SQLDatabase (see
    sage.databases.database.py).  The display functions and get_list create
    their own queries, but it is also possible to query the database by
    constructing either a GenericSQLQuery or a SQLQuery.

    The database contains all unlabeled graphs with 7 or fewer nodes.
    This class will also interface with the optional database package
    containing all unlabeled graphs with 8 or fewer nodes.
    The database(s) consists of five tables, and has the structure given
    by the skeleton:
    \begin{verbatim}
    {'aut_grp': {'aut_grp_size': {'index': True,
                                    'primary_key': False,
                                    'sql': 'INTEGER'},
                  'edge_transitive': {'index': True,
                                       'primary_key': False,
                                       'sql': 'BOOLEAN'},
                  'graph_id': {'index': False,
                                'primary_key': False,
                                'sql': 'INTEGER'},
                  'num_fixed_points': {'index': True,
                                        'primary_key': False,
                                        'sql': 'INTEGER'},
                  'num_orbits': {'index': True,
                                  'primary_key': False,
                                  'sql': 'INTEGER'},
                  'vertex_transitive': {'index': True,
                                         'primary_key': False,
                                         'sql': 'BOOLEAN'}},
     'degrees': {'average_degree': {'index': True,
                                      'primary_key': False,
                                      'sql': 'REAL'},
                  'degree_sequence': {'index': False,
                                       'primary_key': False,
                                       'sql': 'INTEGER'},
                  'degrees_sd': {'index': True,
                                  'primary_key': False,
                                  'sql': 'REAL'},
                  'graph_id': {'index': False,
                                'primary_key': False,
                                'sql': 'INTEGER'},
                  'max_degree': {'index': True,
                                  'primary_key': False,
                                  'sql': 'INTEGER'},
                  'min_degree': {'index': True,
                                  'primary_key': False,
                                  'sql': 'INTEGER'},
                  'regular': {'index': True,
                               'primary_key': False,
                               'sql': 'BOOLEAN'}},
     'graph_data': {'complement_graph6': {'index': True,
                                            'primary_key': False,
                                            'sql': 'TEXT'},
                     'eulerian': {'index': True,
                                   'primary_key': False,
                                   'sql': 'BOOLEAN'},
                     'graph6': {'index': True,
                                 'primary_key': False,
                                 'sql': 'TEXT'},
                     'graph_id': {'index': True,
                                   'primary_key': False,
                                   'sql': 'INTEGER'},
                     'lovasz_number': {'index': True,
                                        'primary_key': False,
                                        'sql': 'REAL'},
                     'num_cycles': {'index': True,
                                     'primary_key': False,
                                     'sql': 'INTEGER'},
                     'num_edges': {'index': True,
                                    'primary_key': False,
                                    'sql': 'INTEGER'},
                     'num_hamiltonian_cycles': {'index': True,
                                                 'primary_key': False,
                                                 'sql': 'INTEGER'},
                     'num_vertices': {'index': True,
                                       'primary_key': False,
                                       'sql': 'INTEGER'},
                     'perfect': {'index': True,
                                  'primary_key': False,
                                  'sql': 'BOOLEAN'},
                     'planar': {'index': True,
                                 'primary_key': False,
                                 'sql': 'BOOLEAN'}},
     'misc': {'clique_number': {'index': True,
                                  'primary_key': False,
                                  'sql': 'INTEGER'},
               'diameter': {'index': True,
                             'primary_key': False,
                             'sql': 'INTEGER'},
               'edge_connectivity': {'index': True,
                                      'primary_key': False,
                                      'sql': 'BOOLEAN'},
               'girth': {'index': True, 'primary_key': False, 'sql': 'INTEGER'},
               'graph_id': {'index': False,
                             'primary_key': False,
                             'sql': 'INTEGER'},
               'independence_number': {'index': True,
                                        'primary_key': False,
                                        'sql': 'INTEGER'},
               'induced_subgraphs': {'index': True,
                                      'primary_key': False,
                                      'sql': 'TEXT'},
               'min_vertex_cover_size': {'index': True,
                                          'primary_key': False,
                                          'sql': 'INTEGER'},
               'num_components': {'index': True,
                                   'primary_key': False,
                                   'sql': 'INTEGER'},
               'num_cut_vertices': {'index': True,
                                     'primary_key': False,
                                     'sql': 'INTEGER'},
               'num_spanning_trees': {'index': True,
                                       'primary_key': False,
                                       'sql': 'INTEGER'},
               'radius': {'index': True,
                           'primary_key': False,
                           'sql': 'INTEGER'},
               'vertex_connectivity': {'index': True,
                                        'primary_key': False,
                                        'sql': 'BOOLEAN'}},
     'spectrum': {'eigenvalues_sd': {'index': True,
                                       'primary_key': False,
                                       'sql': 'REAL'},
                   'energy': {'index': True,
                               'primary_key': False,
                               'sql': 'REAL'},
                   'graph_id': {'index': False,
                                 'primary_key': False,
                                 'sql': 'INTEGER'},
                   'max_eigenvalue': {'index': True,
                                       'primary_key': False,
                                       'sql': 'REAL'},
                   'min_eigenvalue': {'index': True,
                                       'primary_key': False,
                                       'sql': 'REAL'},
                   'spectrum': {'index': False,
                                 'primary_key': False,
                                 'sql': 'TEXT'}}}
    \end{verbatim}

    The automatically generated queries will search graph_data
    automatically, and the other tables will be searched as necessary.
    (Display functions require that all tablesbe searched).

USE:
    The tables are associated by the unique primary key graph_id (int).

    For the query generating functions (display functions and get_list),
    the query parameter allows the user to input any GenericSQLQuery
    associated with this database.  Otherwise, the user can enter
    the individual parameters and the sqlite call will be generated by
    the function.

    The properties currently used as search parameters are:
    \begin{verbatim}
        Table: aut_grp
            - aut_grp_size (The size of the automorphism group - Integer)
            - edge_transitive (Boolean)
            - num_fixed_points (Integer)
            - num_orbits (Integer)
            - vertex_transitive (Boolean)
        Table: degrees
            - average_degree (Real)
            - degree_sequence (Integer)
            - degrees_sd (Standard Deviation of degrees - Real)
            - max_degree (Integer)
            - min_degree (Integer)
            - regular (Boolean)
        Table: graph_data
            - complement_graph6 (graph6 canonical label of the complement
              graph - String)
            - eulerian (Boolean)
            - graph6 (canonical label - String)
            - lovasz_number (Real)
            - num_cycles (Integer)
            - num_edges (Integer)
            - num_hamiltonian_cycles (Integer)
            - num_vertices (Integer)
            - perfect (Boolean)
            - planar (Boolean)
        Table: misc
            - clique_number (Integer)
            - diameter (Real)
            - edge_connectivity (Integer)
            - girth (Integer)
            - independence_number (Integer)
            - induced_subgraphs (canonical label, use regexp - String)
            - min_vertex_cover_size (Integer)
            - num_components (Integer)
            - num_cut_vertices (Integer)
            - num_spanning_trees (Integer)
            - radius (Real)
            - vertex_connectivity (Integer)
        Table: spectrum
            - eigenvalues_sd (Standard Deviation of eigenvalues - Real)
            - energy (Real)
            - max_eigenvalue (Real)
            - min_eigenvalue (Real)
            - spectrum (String)
    \end{verbatim}

VISUALIZATION:

    Beyond the typical show function, there are three options for displaying
    the results of a query.  When running the notebook, each of these functions
    displays an image of the graph and it's (canonical label) graph6 string
    in an html results table.

    \begin{verbatim}
        - display_all (displays all the database properties in the results
          table).
        - display_tables (displays all the properties in the database tables
          that are listed by the user).
        - display_properties (displays all the individual properties
          specified by the user).
    \end{verbatim}

AUTHORS:
    -- Emily A. Kirkman (2007-07-23): inherits GenericSQLDatabase, also added
                                    classes: GraphQuery and GenericGraphQuery
    -- Emily A. Kirkman (2007-05-11): initial sqlite version
    -- Emily A. Kirkman (2007-02-13): initial version (non-sqlite)

REFERENCES:
    -- Data provided by Jason Grout (Brigham Young University).
       [Online] Available: http://math.byu.edu/~grout/graphs/

"""

################################################################################
#           Copyright (C) 2007 Emily A. Kirkman
#
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
################################################################################
import graph
import re
from sqlite3 import dbapi2 as sqlite
import os
from sage.databases.database import GenericSQLDatabase, GenericSQLQuery, SQLQuery

from sage.databases.db import DB_HOME
dblocation = DB_HOME + '/graphs/graphs.db'

def regexp(expr, item):
    """
    Function to define regular expressions in pysqlite.
    Returns 1 if parameter `item` matches the regular expression parameter `expr`.
    Returns 0 otherwise (i.e.: no match).

    REFERENCES:
        Gerhard Haring. [Online] Available: http://lists.initd.org/pipermail/pysqlite/2005-November/000253.html
    """
    r = re.compile(expr)
    return r.match(item) is not None

class GenericGraphQuery(GenericSQLQuery):
    """
    A query for a GraphDatabase.

    INPUT:
        database -- the GraphDatabase instance to query
        query_string -- a string representing the SQL query
        param_tuple -- a tuple of strings - what to replace question marks in
                query_string with (optional, but a good idea)

    NOTE:
        This query class is generally intended for developers and more
        advanced users.  It allows you to execute any query, and so may be
        considered unsafe...

        See GraphDatabase class docstrings or enter:

        sage: G = GraphDatabase()
        sage.: G.get_skeleton()

        to see the underlying structure of the database.  Also see
        GenericSQLQuery in sage.databases.database for more info and a tutorial.

        A piece of advice about '?' and param_tuple:
        It is generally considered safer to query with a '?' in place of
        each value parameter, and using a second argument (a tuple of strings)
        in a call to the sqlite database.  Successful use of the param_tuple
        argument is exemplified:

        sage: G = GraphDatabase()
        sage: q = 'select graph_id,graph6,num_vertices,num_edges from graph_data where graph_id<=(?) and num_vertices=(?)'
        sage: param = (22,5)
        sage: Q = GenericSQLQuery(G,q,param)
        sage: Q.show()
        graph_id             graph6               num_vertices         num_edges
        --------------------------------------------------------------------------------
        18                   D??                  5                    0
        19                   D?C                  5                    1
        20                   D?K                  5                    2
        21                   D@O                  5                    2
        22                   D?[                  5                    3

    """
    def __init__(self, database, query_string, param_tuple=None):
        if not isinstance(database, GraphDatabase):
            raise TypeError('%s is not a valid GraphDatabase'%database)
        GenericSQLQuery.__init__(self,database,query_string,param_tuple)

    def show(self, max_field_size=20, html_table=False, with_picture=False):
        """
        Displays the results of a query in table format.

        INPUT:
            max_field_size -- width of fields in command prompt version
            html_table -- whether or not to draw table in html format
            with_picture -- whether or not to display results with a picture
                of the graph

        EXAMPLES:
        TODO
        """
        if picture:
            from sage.plot.plot import plot

            s = (self.__query_string__).lower()
            s = re.sub('select ','select graph_data.graph6,',s)
            s = re.sub(' from .* where ', ' from graph_data inner join aut_grp on aut_grp.graph_id=graph_data.graph_id inner join degrees on degrees.graph_id=graph_data.graph_id inner join misc on misc.graph_id=graph_data.graph_id inner join spectrum on spectrum.graph_id=graph_data.graph_id where ',s)

            try:
                cur = self.__database__.__connection__.cursor()
                if self.__param_tuple__ is not None:
                    tup = []
                    # make it a tuple of strings:
                    for i in range (len(self.__param_tuple__)):
                        tup.append(str(self.__param_tuple__[i]))
                    cur.execute(s, tuple(tup))
                else:
                    cur.execute(s)
                b = cur.fetchall()
            except:
                raise RuntimeError('Failure to fetch query.')

            # Picture drawn here
            graph6list = []
            for i in range (len(b)):
                graph6 = str(b[i][0])
                g = graph.Graph('%s'%graph6)
                graph6list.append(g.graph6_string())
                p = g.plot(layout='circular', vertex_size=30, vertex_labels=False, graph_border=False)
                p.save('%s.png'%i, figsize=[1,1])

            print '<html><table bgcolor=lightgrey cellpadding=0><tr>'
            print '<td bgcolor=white align=center> Picture </td>'
            for des in cur.description[1:]:
                print '<td bgcolor=white align=center> ' + des[0] + ' </td>'
            print '</tr>'
            field_indices = range(len(cur.description))
            count = 0
            for row in b:
                print '<tr><td bgcolor=white align=center><img src="cell://%d.png"></td>'%count
                for index in field_indices[1:]:
                    print '<td bgcolor=white align=center> ' + str(row[index]) + ' </td>'
                print '</tr>'
                count += 1
            print '</table></html>'

        else:
            GenericSQLQuery.show(self, max_field_size, html_table)

class GraphQuery(SQLQuery):
    """
    A query for an instance of GraphDatabase.  This class nicely wraps
    the SQLQuery class located in sage.databases.database.py to make the
    query constraints intuitive and with as many predefinitions as posible.
    (i.e.: since it has to be a GraphDatabase, we already know the table
    structure and types; and since it is immutable, we can treat these as
    a guarantee).

    NOTE:
        SQLQuery functions are available for GraphQuery.

    TUTORIAL:
    TODO

    INPUT:
        database -- the instance of GraphDatabase to apply query to
        query -- (GenericSQLQuery) A sqlite query for graphs.db (See examples below).
        layout -- (String) The layout option for the graph image.  Options include:
                           'circular' -- plots the graph with vertices evenly
                                         distributed on a circle
                           'spring' -- uses the traditional spring layout
        aut_grp_size -- (Integer) The desired size of the automorphism group.
                        (List) Format: [<String>,<Integer>] WHERE the first
                                       entry represents an inequality:
                                       '=','>','<','>=','<='
        average_degree -- (Real) The desired average degree.
                          (List) Format: [<String>,<Real>] WHERE the first
                                         entry represents an inequality:
                                         '=','>','<','>=','<='
        clique_number -- (Integer) The desired clique number.
                         (List) Format: [<String>,<Integer>] WHERE the first
                                        entry represents an inequality:
                                        '=','>','<','>=','<='
        complement_graph6 -- (String) A graph6 string isomorphic to the
                                      desired complement graph.
                             (List) A list of graph6 strings.  Will search
                                    for graphs with complement isomorphic to
                                    any string in the list.
        degree_sequence -- (Integer) The desired sequence of degrees.
                                     (Ordered highest to lowest).
        degrees_sd -- (Real) The desired standard deviation of degrees.
                      (List) Format: [<String>,<Real>] WHERE the first
                                     entry represents an inequality:
                                     '=','>','<','>=','<='
        diameter -- (Real) The desired diameter.
                    (List) Format: [<String>,<Real>] WHERE the first
                                    entry represents an inequality:
                                    '=','>','<','>=','<='
        edge_connectivity -- (Integer) The desired edge connectivity.
                             (List) Format: [<String>,<Integer>] WHERE the first
                                    entry represents an inequality:
                                    '=','>','<','>=','<='
        edge_transitive -- (Boolean)
        eigenvalues_sd -- (Real) The desired standard deviation of eigenvalues.
                          (List) Format: [<String>,<Real>] WHERE the first
                                 entry represents an inequality:
                                 '=','>','<','>=','<='
        energy -- (Real) The desired energy.
                  (List) Format: [<String>,<Real>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
        eulerian -- (Boolean)
        girth -- (Integer) The desired girth.
                 (List) Format: [<String>,<Integer>] WHERE the first entry
                                represents an inequality: '=','>','<','>=','<='
        graph6 -- (String) A graph6 string isomorphic to the desired graph.
                  (List) A list of graph6 strings.  Will search for graphs
                         isomorphic to any string in the list.
        independence_number -- (Integer) The desired independence number.
                               (List) Format: [<String>,<Integer>] WHERE the
                               first entry represents an inequality:
                               '=','>','<','>=','<='
        induced_subgraphs -- (String) graph6 string isomorphic to desired subgraph.
                             (List) Format options:
                                    1. ['one_of',<String>,...,<String>]
                                       Will search for graphs containing a subgraph
                                       isomorphic to any of the graph6 strings in
                                       the list.
                                    2. ['all_of',<String>,...,<String>]
                                       Will search for graphs containing a subgraph
                                       isomorphic to each of the graph6 strings in
                                       the list.
        lovasz_number -- (Real) The desired lovasz number.
                         (List) Format: [<String>,<Real>] WHERE the first entry
                                represents an inequality: '=','>','<','>=','<='
        max_degree -- (Integer) The desired maximum degree.
                      (List) Format: [<String>,<Integer>] WHERE the first entry
                             represents an inequality: '=','>','<','>=','<='
        max_eigenvalue -- (Real) The desired maximum eigenvalue.
                          (List) Format: [<String>,<Real>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
        min_degree -- (Integer) The desired minimum degree.
                      (List) Format: [<String>,<Integer>] WHERE the first entry
                             represents an inequality: '=','>','<','>=','<='
        min_eigenvalue -- (Real) The desired minimum eigenvalue.
                          (List) Format: [<String>,<Real>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
        min_vertex_cover_size -- (Integer) The desired minimum vertex cover size.
                                 (List) Format: [<String>,<Integer>] WHERE the
                                        first entry represents an inequality:
                                        '=','>','<','>=','<='
        num_components -- (Integer) The desired number of components.
                          (List) Format: [<String>,<Integer>] WHERE the first
                                 entry represents an inequality:
                                 '=','>','<','>=','<='
        num_cut_vertices -- (Integer) The desired number of cut vertices.
                            (List) Format: [<String>,<Integer>] WHERE the first
                                 entry represents an inequality:
                                 '=','>','<','>=','<='
        num_cycles -- (Integer) The desired number of cycles.
                      (List) Format: [<String>,<Integer>] WHERE the first entry
                             represents an inequality: '=','>','<','>=','<='
        num_edges -- (Integer) The desired number of edges.
                     (List) Format: [<String>,<Integer>] WHERE the first entry
                            represents an inequality: '=','>','<','>=','<='
        num_fixed_points -- (Integer) The desired number of fixed points.
                            (List) Format: [<String>,<Integer>] WHERE the first
                                   entry represents an inequality:
                                   '=','>','<','>=','<='
        num_hamiltonian_cycles -- (Integer) The desired number of hamiltonian cycles.
                                  (List) Format: [<String>,<Integer>] WHERE the first
                                         entry represents an inequality:
                                         '=','>','<','>=','<='
        num_orbits -- (Integer) The desired number of orbits.
                      (List) Format: [<String>,<Integer>] WHERE the first entry
                             represents an inequality: '=','>','<','>=','<='
        num_spanning_trees -- (Integer) The desired number of spanning trees.
                              (List) Format: [<String>,<Integer>] WHERE the first
                                     entry represents an inequality:
                                     '=','>','<','>=','<='
        num_vertices -- (Integer) The desired number of vertices.
                        (List) Format: [<String>,<Integer>] WHERE the first entry
                        represents an inequality: '=','>','<','>=','<='
        perfect -- (Boolean)
        planar -- (Boolean)
        radius -- (Integer) The desired radius.
                  (List) Format: [<String>,<Integer>] WHERE the first entry represents
                         an inequality: '=','>','<','>=','<='
        regular -- (Boolean)
        spectrum -- (String) The desired spectrum.  (Ordered highest to lowest,
                             delimited by ', ' and rounded to 6 decimal places).
        vertex_connectivity -- (Integer) The desired vertex connectivity.
                               (List) Format: [<String>,<Integer>] WHERE the first
                                      entry represents an inequality:
                                      '=','>','<','>=','<='
        vertex_transitive -- (Boolean)
    """
    def __init__(self, database, graph6=None, num_vertices=None, \
                    num_edges=None, num_cycles=None, num_hamiltonian_cycles=None, \
                    eulerian=None, planar=None, perfect=None, lovasz_number=None, \
                    complement_graph6=None, aut_grp_size=None, num_orbits=None, \
                    num_fixed_points=None, vertex_transitive=None, edge_transitive=None, \
                    degree_sequence=None, min_degree=None, max_degree=None, \
                    average_degree=None, degrees_sd=None, regular=None, \
                    vertex_connectivity=None, edge_connectivity=None, num_components=None, \
                    girth=None, radius=None, diameter=None, clique_number=None, \
                    independence_number=None, num_cut_vertices=None, \
                    min_vertex_cover_size=None, num_spanning_trees=None, \
                    induced_subgraphs=None, spectrum=None, min_eigenvalue=None, \
                    max_eigenvalue=None, eigenvalues_sd=None, energy=None):

        SQLQuery.__init__(self,database)

        query = __query_string__(graph6=graph6, num_vertices=num_vertices, num_edges=num_edges, num_cycles=num_cycles, num_hamiltonian_cycles=num_hamiltonian_cycles, eulerian=eulerian, planar=planar, perfect=perfect, lovasz_number=lovasz_number, complement_graph6=complement_graph6, aut_grp_size=aut_grp_size, num_orbits=num_orbits, num_fixed_points=num_fixed_points, vertex_transitive=vertex_transitive, edge_transitive=edge_transitive, degree_sequence=degree_sequence, min_degree=min_degree, max_degree=max_degree, average_degree=average_degree, degrees_sd=degrees_sd, regular=regular, vertex_connectivity=vertex_connectivity, edge_connectivity=edge_connectivity, num_components=num_components, girth=girth, radius=radius, diameter=diameter, clique_number=clique_number, independence_number=independence_number, num_cut_vertices=num_cut_vertices, min_vertex_cover_size=min_vertex_cover_size, num_spanning_trees=num_spanning_trees, induced_subgraphs=induced_subgraphs, spectrum=spectrum, min_eigenvalue=min_eigenvalue, max_eigenvalue=max_eigenvalue, eigenvalues_sd=eigenvalues_sd, energy=energy)
        query = re.sub('INNER JOIN .* WHERE','INNER JOIN aut_grp on aut_grp.graph_id=graph_data.graph_id INNER JOIN degrees on degrees.graph_id=graph_data.graph_id INNER JOIN misc on misc.graph_id=graph_data.graph_id INNER JOIN spectrum on spectrum.graph_id=graph_data.graph_id WHERE',query)
        query = re.sub('FROM graph_data WHERE','FROM graph_data INNER JOIN aut_grp on aut_grp.graph_id=graph_data.graph_id INNER JOIN degrees on degrees.graph_id=graph_data.graph_id INNER JOIN misc on misc.graph_id=graph_data.graph_id INNER JOIN spectrum on spectrum.graph_id=graph_data.graph_id WHERE',query)
        query = re.sub('SELECT graph_data.graph6','SELECT graph_data.graph6,graph_data.num_vertices,degrees.regular,aut_grp.aut_grp_size,graph_data.num_edges,degrees.min_degree,aut_grp.num_orbits,graph_data.num_cycles,degrees.max_degree,aut_grp.num_fixed_points,graph_data.num_hamiltonian_cycles,degrees.average_degree,aut_grp.vertex_transitive,graph_data.eulerian,degrees.degrees_sd,aut_grp.edge_transitive,graph_data.planar,degrees.degree_sequence,misc.vertex_connectivity,graph_data.perfect,spectrum.min_eigenvalue,misc.edge_connectivity,graph_data.lovasz_number,misc.girth,spectrum.max_eigenvalue,misc.num_cut_vertices,misc.independence_number,misc.radius,spectrum.eigenvalues_sd,misc.min_vertex_cover_size,misc.clique_number,misc.diameter,spectrum.energy,misc.num_spanning_trees,misc.num_components,graph_data.complement_graph6,spectrum.spectrum,misc.induced_subgraphs',query)

        self.__query_string__ = query

class GraphDatabase(GenericSQLDatabase):
    r"""
    Graph Database

    INFO:

        This class interfaces with the sqlite database graphs.db.  It is an
        immutable database that inherits from SQLDatabase (see
        sage.databases.database.py).  The display functions and get_list create
        their own queries, but it is also possible to query the database by
        constructing either a GenericSQLQuery or a SQLQuery.

        The database contains all unlabeled graphs with 7 or fewer nodes.
        This class will also interface with the optional database package
        containing all unlabeled graphs with 8 or fewer nodes.
        The database(s) consists of five tables, and has the structure given
        by the skeleton:
        \begin{verbatim}
        {'aut_grp': {'aut_grp_size': {'index': True,
                                        'primary_key': False,
                                        'sql': 'INTEGER'},
                      'edge_transitive': {'index': True,
                                           'primary_key': False,
                                           'sql': 'BOOLEAN'},
                      'graph_id': {'index': False,
                                    'primary_key': False,
                                    'sql': 'INTEGER'},
                      'num_fixed_points': {'index': True,
                                            'primary_key': False,
                                            'sql': 'INTEGER'},
                      'num_orbits': {'index': True,
                                      'primary_key': False,
                                      'sql': 'INTEGER'},
                      'vertex_transitive': {'index': True,
                                             'primary_key': False,
                                             'sql': 'BOOLEAN'}},
         'degrees': {'average_degree': {'index': True,
                                          'primary_key': False,
                                          'sql': 'REAL'},
                      'degree_sequence': {'index': False,
                                           'primary_key': False,
                                           'sql': 'INTEGER'},
                      'degrees_sd': {'index': True,
                                      'primary_key': False,
                                      'sql': 'REAL'},
                      'graph_id': {'index': False,
                                    'primary_key': False,
                                    'sql': 'INTEGER'},
                      'max_degree': {'index': True,
                                      'primary_key': False,
                                      'sql': 'INTEGER'},
                      'min_degree': {'index': True,
                                      'primary_key': False,
                                      'sql': 'INTEGER'},
                      'regular': {'index': True,
                                   'primary_key': False,
                                   'sql': 'BOOLEAN'}},
         'graph_data': {'complement_graph6': {'index': True,
                                                'primary_key': False,
                                                'sql': 'TEXT'},
                         'eulerian': {'index': True,
                                       'primary_key': False,
                                       'sql': 'BOOLEAN'},
                         'graph6': {'index': True,
                                     'primary_key': False,
                                     'sql': 'TEXT'},
                         'graph_id': {'index': True,
                                       'primary_key': False,
                                       'sql': 'INTEGER'},
                         'lovasz_number': {'index': True,
                                            'primary_key': False,
                                            'sql': 'REAL'},
                         'num_cycles': {'index': True,
                                         'primary_key': False,
                                         'sql': 'INTEGER'},
                         'num_edges': {'index': True,
                                        'primary_key': False,
                                        'sql': 'INTEGER'},
                         'num_hamiltonian_cycles': {'index': True,
                                                     'primary_key': False,
                                                     'sql': 'INTEGER'},
                         'num_vertices': {'index': True,
                                           'primary_key': False,
                                           'sql': 'INTEGER'},
                         'perfect': {'index': True,
                                      'primary_key': False,
                                      'sql': 'BOOLEAN'},
                         'planar': {'index': True,
                                     'primary_key': False,
                                     'sql': 'BOOLEAN'}},
         'misc': {'clique_number': {'index': True,
                                      'primary_key': False,
                                      'sql': 'INTEGER'},
                   'diameter': {'index': True,
                                 'primary_key': False,
                                 'sql': 'INTEGER'},
                   'edge_connectivity': {'index': True,
                                          'primary_key': False,
                                          'sql': 'BOOLEAN'},
                   'girth': {'index': True, 'primary_key': False, 'sql': 'INTEGER'},
                   'graph_id': {'index': False,
                                 'primary_key': False,
                                 'sql': 'INTEGER'},
                   'independence_number': {'index': True,
                                            'primary_key': False,
                                            'sql': 'INTEGER'},
                   'induced_subgraphs': {'index': True,
                                          'primary_key': False,
                                          'sql': 'TEXT'},
                   'min_vertex_cover_size': {'index': True,
                                              'primary_key': False,
                                              'sql': 'INTEGER'},
                   'num_components': {'index': True,
                                       'primary_key': False,
                                       'sql': 'INTEGER'},
                   'num_cut_vertices': {'index': True,
                                         'primary_key': False,
                                         'sql': 'INTEGER'},
                   'num_spanning_trees': {'index': True,
                                           'primary_key': False,
                                           'sql': 'INTEGER'},
                   'radius': {'index': True,
                               'primary_key': False,
                               'sql': 'INTEGER'},
                   'vertex_connectivity': {'index': True,
                                            'primary_key': False,
                                            'sql': 'BOOLEAN'}},
         'spectrum': {'eigenvalues_sd': {'index': True,
                                           'primary_key': False,
                                           'sql': 'REAL'},
                       'energy': {'index': True,
                                   'primary_key': False,
                                   'sql': 'REAL'},
                       'graph_id': {'index': False,
                                     'primary_key': False,
                                     'sql': 'INTEGER'},
                       'max_eigenvalue': {'index': True,
                                           'primary_key': False,
                                           'sql': 'REAL'},
                       'min_eigenvalue': {'index': True,
                                           'primary_key': False,
                                           'sql': 'REAL'},
                       'spectrum': {'index': False,
                                     'primary_key': False,
                                     'sql': 'TEXT'}}}
        \end{verbatim}

        The automatically generated queries will search graph_data
        automatically, and the other tables will be searched as necessary.
        (Display functions require that all tablesbe searched).

    USE:
        The tables are associated by the unique primary key graph_id (int).

        For the query generating functions (display functions and get_list),
        the query parameter allows the user to input any GenericSQLQuery
        associated with this database.  Otherwise, the user can enter
        the individual parameters and the sqlite call will be generated by
        the function.

        The properties currently used as search parameters are:
        \begin{verbatim}
            Table: aut_grp
                - aut_grp_size (The size of the automorphism group - Integer)
                - edge_transitive (Boolean)
                - num_fixed_points (Integer)
                - num_orbits (Integer)
                - vertex_transitive (Boolean)
            Table: degrees
                - average_degree (Real)
                - degree_sequence (Integer)
                - degrees_sd (Standard Deviation of degrees - Real)
                - max_degree (Integer)
                - min_degree (Integer)
                - regular (Boolean)
            Table: graph_data
                - complement_graph6 (graph6 canonical label of the complement
                  graph - String)
                - eulerian (Boolean)
                - graph6 (canonical label - String)
                - lovasz_number (Real)
                - num_cycles (Integer)
                - num_edges (Integer)
                - num_hamiltonian_cycles (Integer)
                - num_vertices (Integer)
                - perfect (Boolean)
                - planar (Boolean)
            Table: misc
                - clique_number (Integer)
                - diameter (Real)
                - edge_connectivity (Integer)
                - girth (Integer)
                - independence_number (Integer)
                - induced_subgraphs (canonical label, use regexp - String)
                - min_vertex_cover_size (Integer)
                - num_components (Integer)
                - num_cut_vertices (Integer)
                - num_spanning_trees (Integer)
                - radius (Real)
                - vertex_connectivity (Integer)
            Table: spectrum
                - eigenvalues_sd (Standard Deviation of eigenvalues - Real)
                - energy (Real)
                - max_eigenvalue (Real)
                - min_eigenvalue (Real)
                - spectrum (String)
        \end{verbatim}

    VISUALIZATION:

        Beyond the typical show function, there are three options for displaying
        the results of a query.  When running the notebook, each of these functions
        displays an image of the graph and it's (canonical label) graph6 string
        in an html results table.

        \begin{verbatim}
            - display_all (displays all the database properties in the results
              table).
            - display_tables (displays all the properties in the database tables
              that are listed by the user).
            - display_properties (displays all the individual properties
              specified by the user).
        \end{verbatim}

    REFERENCES:
        -- Data provided by Jason Grout (Brigham Young University).
           [Online] Available: http://math.byu.edu/~grout/graphs/

    """

    def __init__(self):
        GenericSQLDatabase.__init__(self,dblocation)

    def display_all(self, query=None, layout='circular', graph6=None, num_vertices=None, \
                    num_edges=None, num_cycles=None, num_hamiltonian_cycles=None, \
                    eulerian=None, planar=None, perfect=None, lovasz_number=None, \
                    complement_graph6=None, aut_grp_size=None, num_orbits=None, \
                    num_fixed_points=None, vertex_transitive=None, edge_transitive=None, \
                    degree_sequence=None, min_degree=None, max_degree=None, \
                    average_degree=None, degrees_sd=None, regular=None, \
                    vertex_connectivity=None, edge_connectivity=None, num_components=None, \
                    girth=None, radius=None, diameter=None, clique_number=None, \
                    independence_number=None, num_cut_vertices=None, \
                    min_vertex_cover_size=None, num_spanning_trees=None, \
                    induced_subgraphs=None, spectrum=None, min_eigenvalue=None, \
                    max_eigenvalue=None, eigenvalues_sd=None, energy=None):
        r"""
        Displays the results of a query in a table, including all stored
        properties and an image for each graph.

        INPUT:
            query -- (GenericSQLQuery) A sqlite query for graphs.db (See examples below).
            layout -- (String) The layout option for the graph image.  Options include:
                               'circular' -- plots the graph with vertices evenly
                                             distributed on a circle
                               'spring' -- uses the traditional spring layout
            aut_grp_size -- (Integer) The desired size of the automorphism group.
                            (List) Format: [<String>,<Integer>] WHERE the first
                                           entry represents an inequality:
                                           '=','>','<','>=','<='
            average_degree -- (Real) The desired average degree.
                              (List) Format: [<String>,<Real>] WHERE the first
                                             entry represents an inequality:
                                             '=','>','<','>=','<='
            clique_number -- (Integer) The desired clique number.
                             (List) Format: [<String>,<Integer>] WHERE the first
                                            entry represents an inequality:
                                            '=','>','<','>=','<='
            complement_graph6 -- (String) A graph6 string isomorphic to the
                                          desired complement graph.
                                 (List) A list of graph6 strings.  Will search
                                        for graphs with complement isomorphic to
                                        any string in the list.
            degree_sequence -- (Integer) The desired sequence of degrees.
                                         (Ordered highest to lowest).
            degrees_sd -- (Real) The desired standard deviation of degrees.
                          (List) Format: [<String>,<Real>] WHERE the first
                                         entry represents an inequality:
                                         '=','>','<','>=','<='
            diameter -- (Real) The desired diameter.
                        (List) Format: [<String>,<Real>] WHERE the first
                                        entry represents an inequality:
                                        '=','>','<','>=','<='
            edge_connectivity -- (Integer) The desired edge connectivity.
                                 (List) Format: [<String>,<Integer>] WHERE the first
                                        entry represents an inequality:
                                        '=','>','<','>=','<='
            edge_transitive -- (Boolean)
            eigenvalues_sd -- (Real) The desired standard deviation of eigenvalues.
                              (List) Format: [<String>,<Real>] WHERE the first
                                     entry represents an inequality:
                                     '=','>','<','>=','<='
            energy -- (Real) The desired energy.
                      (List) Format: [<String>,<Real>] WHERE the first entry
                                     represents an inequality: '=','>','<','>=','<='
            eulerian -- (Boolean)
            girth -- (Integer) The desired girth.
                     (List) Format: [<String>,<Integer>] WHERE the first entry
                                    represents an inequality: '=','>','<','>=','<='
            graph6 -- (String) A graph6 string isomorphic to the desired graph.
                      (List) A list of graph6 strings.  Will search for graphs
                             isomorphic to any string in the list.
            independence_number -- (Integer) The desired independence number.
                                   (List) Format: [<String>,<Integer>] WHERE the
                                   first entry represents an inequality:
                                   '=','>','<','>=','<='
            induced_subgraphs -- (String) graph6 string isomorphic to desired subgraph.
                                 (List) Format options:
                                        1. ['one_of',<String>,...,<String>]
                                           Will search for graphs containing a subgraph
                                           isomorphic to any of the graph6 strings in
                                           the list.
                                        2. ['all_of',<String>,...,<String>]
                                           Will search for graphs containing a subgraph
                                           isomorphic to each of the graph6 strings in
                                           the list.
            lovasz_number -- (Real) The desired lovasz number.
                             (List) Format: [<String>,<Real>] WHERE the first entry
                                    represents an inequality: '=','>','<','>=','<='
            max_degree -- (Integer) The desired maximum degree.
                          (List) Format: [<String>,<Integer>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
            max_eigenvalue -- (Real) The desired maximum eigenvalue.
                              (List) Format: [<String>,<Real>] WHERE the first entry
                                     represents an inequality: '=','>','<','>=','<='
            min_degree -- (Integer) The desired minimum degree.
                          (List) Format: [<String>,<Integer>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
            min_eigenvalue -- (Real) The desired minimum eigenvalue.
                              (List) Format: [<String>,<Real>] WHERE the first entry
                                     represents an inequality: '=','>','<','>=','<='
            min_vertex_cover_size -- (Integer) The desired minimum vertex cover size.
                                     (List) Format: [<String>,<Integer>] WHERE the
                                            first entry represents an inequality:
                                            '=','>','<','>=','<='
            num_components -- (Integer) The desired number of components.
                              (List) Format: [<String>,<Integer>] WHERE the first
                                     entry represents an inequality:
                                     '=','>','<','>=','<='
            num_cut_vertices -- (Integer) The desired number of cut vertices.
                                (List) Format: [<String>,<Integer>] WHERE the first
                                     entry represents an inequality:
                                     '=','>','<','>=','<='
            num_cycles -- (Integer) The desired number of cycles.
                          (List) Format: [<String>,<Integer>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
            num_edges -- (Integer) The desired number of edges.
                         (List) Format: [<String>,<Integer>] WHERE the first entry
                                represents an inequality: '=','>','<','>=','<='
            num_fixed_points -- (Integer) The desired number of fixed points.
                                (List) Format: [<String>,<Integer>] WHERE the first
                                       entry represents an inequality:
                                       '=','>','<','>=','<='
            num_hamiltonian_cycles -- (Integer) The desired number of hamiltonian cycles.
                                      (List) Format: [<String>,<Integer>] WHERE the first
                                             entry represents an inequality:
                                             '=','>','<','>=','<='
            num_orbits -- (Integer) The desired number of orbits.
                          (List) Format: [<String>,<Integer>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
            num_spanning_trees -- (Integer) The desired number of spanning trees.
                                  (List) Format: [<String>,<Integer>] WHERE the first
                                         entry represents an inequality:
                                         '=','>','<','>=','<='
            num_vertices -- (Integer) The desired number of vertices.
                            (List) Format: [<String>,<Integer>] WHERE the first entry
                            represents an inequality: '=','>','<','>=','<='
            perfect -- (Boolean)
            planar -- (Boolean)
            radius -- (Integer) The desired radius.
                      (List) Format: [<String>,<Integer>] WHERE the first entry represents
                             an inequality: '=','>','<','>=','<='
            regular -- (Boolean)
            spectrum -- (String) The desired spectrum.  (Ordered highest to lowest,
                                 delimited by ', ' and rounded to 6 decimal places).
            vertex_connectivity -- (Integer) The desired vertex connectivity.
                                   (List) Format: [<String>,<Integer>] WHERE the first
                                          entry represents an inequality:
                                          '=','>','<','>=','<='
            vertex_transitive -- (Boolean)

        EXAMPLES:
        TODO
        The basics:
            sage.: graphs_query.display_all(num_vertices=5,lovasz_number=3.0,\
            ...                             girth=4,radius=2,diameter=3)
            sage.: graphs_query.display_all(layout='spring',num_hamiltonian_cycles=2,\
            ...                             regular=True,perfect=False)
            sage.: graphs_query.display_all(layout='spring',degree_sequence=433211)

        Compare results:
            sage: (Graph('EAMw')).is_isomorphic(Graph('E@NW'))
            False

        Using Inequalities:
            sage.: graphs_query.display_all(layout='circular', min_eigenvalue=['=',-1], \
            ...                             eigenvalues_sd=['<=',1], energy=['>',5])

        The query string:
            sage.: graphs_query.display_all(layout='spring', \
            ...             query='SELECT graph_data.graph6 \
            ...             FROM graph_data WHERE num_vertices<=4 \
            ...             and num_edges>3')
            sage.: graphs_query.display_all(query='SELECT graph_data.graph6 FROM graph_data \
            ...             INNER JOIN degrees on graph_data.graph_id=degrees.graph_id \
            ...             WHERE num_vertices>6 and eulerian=1 and regular=0 and planar=1 \
            ...             and num_cycles<=2')
            sage.: graphs_query.display_all(query="SELECT graph_data.graph6 \
            ...             FROM graph_data INNER JOIN misc on \
            ...             misc.graph_id=graph_data.graph_id WHERE \
            ...             misc.induced_subgraphs regexp '.*E~~w.*'")
        """
        from sage.plot.plot import plot

        if ( query is None):
            param = None
            query = __query_string__(graph6=graph6, num_vertices=num_vertices, num_edges=num_edges, num_cycles=num_cycles, num_hamiltonian_cycles=num_hamiltonian_cycles, eulerian=eulerian, planar=planar, perfect=perfect, lovasz_number=lovasz_number, complement_graph6=complement_graph6, aut_grp_size=aut_grp_size, num_orbits=num_orbits, num_fixed_points=num_fixed_points, vertex_transitive=vertex_transitive, edge_transitive=edge_transitive, degree_sequence=degree_sequence, min_degree=min_degree, max_degree=max_degree, average_degree=average_degree, degrees_sd=degrees_sd, regular=regular, vertex_connectivity=vertex_connectivity, edge_connectivity=edge_connectivity, num_components=num_components, girth=girth, radius=radius, diameter=diameter, clique_number=clique_number, independence_number=independence_number, num_cut_vertices=num_cut_vertices, min_vertex_cover_size=min_vertex_cover_size, num_spanning_trees=num_spanning_trees, induced_subgraphs=induced_subgraphs, spectrum=spectrum, min_eigenvalue=min_eigenvalue, max_eigenvalue=max_eigenvalue, eigenvalues_sd=eigenvalues_sd, energy=energy)
            query = re.sub('INNER JOIN .* WHERE','INNER JOIN aut_grp on aut_grp.graph_id=graph_data.graph_id INNER JOIN degrees on degrees.graph_id=graph_data.graph_id INNER JOIN misc on misc.graph_id=graph_data.graph_id INNER JOIN spectrum on spectrum.graph_id=graph_data.graph_id WHERE',query)
            query = re.sub('FROM graph_data WHERE','FROM graph_data INNER JOIN aut_grp on aut_grp.graph_id=graph_data.graph_id INNER JOIN degrees on degrees.graph_id=graph_data.graph_id INNER JOIN misc on misc.graph_id=graph_data.graph_id INNER JOIN spectrum on spectrum.graph_id=graph_data.graph_id WHERE',query)
            query = re.sub('SELECT graph_data.graph6','SELECT graph_data.graph6,graph_data.num_vertices,degrees.regular,aut_grp.aut_grp_size,graph_data.num_edges,degrees.min_degree,aut_grp.num_orbits,graph_data.num_cycles,degrees.max_degree,aut_grp.num_fixed_points,graph_data.num_hamiltonian_cycles,degrees.average_degree,aut_grp.vertex_transitive,graph_data.eulerian,degrees.degrees_sd,aut_grp.edge_transitive,graph_data.planar,degrees.degree_sequence,misc.vertex_connectivity,graph_data.perfect,spectrum.min_eigenvalue,misc.edge_connectivity,graph_data.lovasz_number,misc.girth,spectrum.max_eigenvalue,misc.num_cut_vertices,misc.independence_number,misc.radius,spectrum.eigenvalues_sd,misc.min_vertex_cover_size,misc.clique_number,misc.diameter,spectrum.energy,misc.num_spanning_trees,misc.num_components,graph_data.complement_graph6,spectrum.spectrum,misc.induced_subgraphs',query)
        else:
            # Deal only with the string:
            param = query.__param_tuple__
            query = query.__query_string__

        cur = (self.__connection__).cursor()
        if param is None:
            a = cur.execute(query)
            b = a.fetchall()
        else:
            tup = []
            # make it a tuple of strings:
            for i in range (len(param)):
                tup.append(str(param[i]))
            exe = cur.execute(query, tuple(tup))
            b = exe.fetchall()

        graph6list = []
        for i in range (len(b)):
            graph6 = str(b[i][0])
            g = graph.Graph('%s'%graph6)
            # The following line is time consuming and should not stay:
            graph6list.append(g.graph6_string())
            p = g.plot(layout=layout, vertex_size=30, vertex_labels=False, graph_border=False)
            p.save('%s.png'%i, figsize=[1,1])

        print "<html>"
        print '<table bgcolor=lightgrey cellpadding=0>'

        from sage.misc.multireplace import multiple_replace
        to_bool = {'0':"False", '1':"True"}

        for i in range(len(b)):
            eul = multiple_replace(to_bool,'%s'%b[i][13])
            reg = multiple_replace(to_bool,'%s'%b[i][2])
            plan = multiple_replace(to_bool,'%s'%b[i][16])
            perf = multiple_replace(to_bool,'%s'%b[i][19])
            vtran = multiple_replace(to_bool,'%s'%b[i][12])
            etran = multiple_replace(to_bool,'%s'%b[i][15])

            print '<tr><td bgcolor=white align=center rowspan="7"><img src="cell://%s.png"><br>%s</td>'%(i,graph6list[i])
            print '<td bgcolor=white align=left><font color=black> Vertices: %s </font></td>'%(b[i][1])
            print '<td bgcolor=white align=left><font color=black> Regular: %s </font></td>'%reg
            print '<td bgcolor=white align=left><font color=black> Aut Group Size: %s </font></td>'%(b[i][3])

            print '</tr><tr><td bgcolor=white align=left><font color=black> Edges: %s </font></td>'%(b[i][4])
            print '<td bgcolor=white align=left><font color=black> Min Degree: %s </font></td>'%(b[i][5])
            print '<td bgcolor=white align=left><font color=black> Orbits: %s </font></td>'%(b[i][6])

            print '</tr><tr><td bgcolor=white align=left><font color=black> Cycles: %s </font></td>'%(b[i][7])
            print '<td bgcolor=white align=left><font color=black> Max Degree: %s </font></td>'%(b[i][8])
            print '<td bgcolor=white align=left><font color=black> Fixed Points: %s </font></td>'%(b[i][9])

            print '</tr><tr><td bgcolor=white align=left><font color=black> Hamiltonian Cycles: %s </font></td>'%(b[i][10])
            print '<td bgcolor=white align=left><font color=black> Average Degree: %s </font></td>'%(b[i][11])
            print '<td bgcolor=white align=left><font color=black> Vertex Transitive: %s </font></td>'%vtran

            print '</tr><tr><td bgcolor=white align=left><font color=black> Eulerian: %s </font></td>'%eul
            print '<td bgcolor=white align=left><font color=black> Degree SD: %s </font></td>'%(b[i][14])
            print '<td bgcolor=white align=left><font color=black> Edge Transitive: %s </font></td>'%etran

            print '</tr><tr><td bgcolor=white align=left><font color=black> Planar: %s </font></td>'%plan
            print '<td bgcolor=white align=left><font color=black> Degree Sequence: %s </font></td>'%(b[i][17])
            print '<td bgcolor=white align=left><font color=black> Vertex Connectivity: %s </font></td>'%(b[i][18])

            print '</tr><tr><td bgcolor=white align=left><font color=black> Perfect: %s </font></td>'%perf
            print '<td bgcolor=white align=left><font color=black> Min Eigenvalue: %s </font></td>'%(b[i][20])
            print '<td bgcolor=white align=left><font color=black> Edge Connectivity: %s </font></td>'%(b[i][21])

            print '</tr><tr><td bgcolor=white align=left><font color=black> Lovasz Number: %s </font></td>'%(b[i][22])
            print '<td bgcolor=white align=left><font color=black> Girth: %s </font></td>'%(b[i][23])
            print '<td bgcolor=white align=left><font color=black> Max Eigenvalue: %s </font></td>'%(b[i][24])
            print '<td bgcolor=white align=left><font color=black> Cut Vertices: %s </font></td>'%(b[i][25])

            print '</tr><tr><td bgcolor=white align=left><font color=black> Independence Number: %s </font></td>'%(b[i][26])
            print '<td bgcolor=white align=left><font color=black> Radius: %s </font></td>'%(b[i][27])
            print '<td bgcolor=white align=left><font color=black> Eigenvalues SD: %s </font></td>'%(b[i][28])
            print '<td bgcolor=white align=left><font color=black> Min Vertex Cover Size: %s </font></td>'%(b[i][29])

            print '</tr><tr><td bgcolor=white align=left><font color=black> Clique Number: %s </font></td>'%(b[i][30])
            print '<td bgcolor=white align=left><font color=black> Diameter: %s </font></td>'%(b[i][31])
            print '<td bgcolor=white align=left><font color=black> Energy: %s </font></td>'%(b[i][32])
            print '<td bgcolor=white align=left><font color=black> Spanning Trees: %s </font></td>'%(b[i][33])

            print '</tr><tr><td bgcolor=white align=left><font color=black> Components: %s </font></td>'%(b[i][34])
            print '<td bgcolor=white align=left><font color=black> Complement: %s </font></td>'%(b[i][35])
            print '<td bgcolor=white align=left colspan="2"><font color=black> Spectrum: %s </font></td>'%(b[i][36])

            print '</tr><tr><td bgcolor=white align=left colspan="4"><font color=black> Induced Subgraphs: %s </font></td></tr>'%(b[i][37])

            if ( i != len(b)-1 ): print '<tr><td bgcolor=lightblue colspan="4" height="3"></td></tr>'

        print "</table></html>"
        return

    def display_tables(self, tables=None, layout='circular', query=None, graph6=None, \
                        num_vertices=None, num_edges=None, num_cycles=None, \
                        num_hamiltonian_cycles=None, eulerian=None, planar=None, perfect=None, \
                        lovasz_number=None, complement_graph6=None, aut_grp_size=None, \
                        num_orbits=None, num_fixed_points=None, vertex_transitive=None, \
                        edge_transitive=None, degree_sequence=None, min_degree=None, \
                        max_degree=None, average_degree=None, degrees_sd=None, regular=None, \
                        vertex_connectivity=None, edge_connectivity=None, num_components=None, \
                        girth=None, radius=None, diameter=None, clique_number=None, \
                        independence_number=None, num_cut_vertices=None, \
                        min_vertex_cover_size=None, num_spanning_trees=None, \
                        induced_subgraphs=None, spectrum=None, min_eigenvalue=None, \
                        max_eigenvalue=None, eigenvalues_sd=None, energy=None):
        r"""
        Displays the results of a query in a table, including all stored
        properties FROM the specified database tables and an image for each graph.

        INPUT:
            query -- (GenericSQLQuery) A sqlite query for graphs.db (See examples below).
            tables -- (List) A list of strings with the exact name (as the
                             database tables) of the tables of properties to
                             display with the results.  Database table names are:
                             'aut_grp','degrees','graph_data','misc','spectrum'
            layout -- (String) The layout option for the graph image.  Options include:
                               'circular' -- plots the graph with vertices evenly
                                             distributed on a circle
                               'spring' -- uses the traditional spring layout
            aut_grp_size -- (Integer) The desired size of the automorphism group.
                            (List) Format: [<String>,<Integer>] WHERE the first
                                           entry represents an inequality:
                                           '=','>','<','>=','<='
            average_degree -- (Real) The desired average degree.
                              (List) Format: [<String>,<Real>] WHERE the first
                                             entry represents an inequality:
                                             '=','>','<','>=','<='
            clique_number -- (Integer) The desired clique number.
                             (List) Format: [<String>,<Integer>] WHERE the first
                                            entry represents an inequality:
                                            '=','>','<','>=','<='
            complement_graph6 -- (String) A graph6 string isomorphic to the
                                          desired complement graph.
                                 (List) A list of graph6 strings.  Will search
                                        for graphs with complement isomorphic to
                                        any string in the list.
            degree_sequence -- (Integer) The desired sequence of degrees.
                                         (Ordered highest to lowest).
            degrees_sd -- (Real) The desired standard deviation of degrees.
                          (List) Format: [<String>,<Real>] WHERE the first
                                         entry represents an inequality:
                                         '=','>','<','>=','<='
            diameter -- (Real) The desired diameter.
                        (List) Format: [<String>,<Real>] WHERE the first
                                        entry represents an inequality:
                                        '=','>','<','>=','<='
            edge_connectivity -- (Integer) The desired edge connectivity.
                                 (List) Format: [<String>,<Integer>] WHERE the first
                                        entry represents an inequality:
                                        '=','>','<','>=','<='
            edge_transitive -- (Boolean)
            eigenvalues_sd -- (Real) The desired standard deviation of eigenvalues.
                              (List) Format: [<String>,<Real>] WHERE the first
                                     entry represents an inequality:
                                     '=','>','<','>=','<='
            energy -- (Real) The desired energy.
                      (List) Format: [<String>,<Real>] WHERE the first entry
                                     represents an inequality: '=','>','<','>=','<='
            eulerian -- (Boolean)
            girth -- (Integer) The desired girth.
                     (List) Format: [<String>,<Integer>] WHERE the first entry
                                    represents an inequality: '=','>','<','>=','<='
            graph6 -- (String) A graph6 string isomorphic to the desired graph.
                      (List) A list of graph6 strings.  Will search for graphs
                             isomorphic to any string in the list.
            independence_number -- (Integer) The desired independence number.
                                   (List) Format: [<String>,<Integer>] WHERE the
                                   first entry represents an inequality:
                                   '=','>','<','>=','<='
            induced_subgraphs -- (String) graph6 string isomorphic to desired subgraph.
                                 (List) Format options:
                                        1. ['one_of',<String>,...,<String>]
                                           Will search for graphs containing a subgraph
                                           isomorphic to any of the graph6 strings in
                                           the list.
                                        2. ['all_of',<String>,...,<String>]
                                           Will search for graphs containing a subgraph
                                           isomorphic to each of the graph6 strings in
                                           the list.
            lovasz_number -- (Real) The desired lovasz number.
                             (List) Format: [<String>,<Real>] WHERE the first entry
                                    represents an inequality: '=','>','<','>=','<='
            max_degree -- (Integer) The desired maximum degree.
                          (List) Format: [<String>,<Integer>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
            max_eigenvalue -- (Real) The desired maximum eigenvalue.
                              (List) Format: [<String>,<Real>] WHERE the first entry
                                     represents an inequality: '=','>','<','>=','<='
            min_degree -- (Integer) The desired minimum degree.
                          (List) Format: [<String>,<Integer>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
            min_eigenvalue -- (Real) The desired minimum eigenvalue.
                              (List) Format: [<String>,<Real>] WHERE the first entry
                                     represents an inequality: '=','>','<','>=','<='
            min_vertex_cover_size -- (Integer) The desired minimum vertex cover size.
                                     (List) Format: [<String>,<Integer>] WHERE the
                                            first entry represents an inequality:
                                            '=','>','<','>=','<='
            num_components -- (Integer) The desired number of components.
                              (List) Format: [<String>,<Integer>] WHERE the first
                                     entry represents an inequality:
                                     '=','>','<','>=','<='
            num_cut_vertices -- (Integer) The desired number of cut vertices.
                                (List) Format: [<String>,<Integer>] WHERE the first
                                     entry represents an inequality:
                                     '=','>','<','>=','<='
            num_cycles -- (Integer) The desired number of cycles.
                          (List) Format: [<String>,<Integer>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
            num_edges -- (Integer) The desired number of edges.
                         (List) Format: [<String>,<Integer>] WHERE the first entry
                                represents an inequality: '=','>','<','>=','<='
            num_fixed_points -- (Integer) The desired number of fixed points.
                                (List) Format: [<String>,<Integer>] WHERE the first
                                       entry represents an inequality:
                                       '=','>','<','>=','<='
            num_hamiltonian_cycles -- (Integer) The desired number of hamiltonian cycles.
                                      (List) Format: [<String>,<Integer>] WHERE the first
                                             entry represents an inequality:
                                             '=','>','<','>=','<='
            num_orbits -- (Integer) The desired number of orbits.
                          (List) Format: [<String>,<Integer>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
            num_spanning_trees -- (Integer) The desired number of spanning trees.
                                  (List) Format: [<String>,<Integer>] WHERE the first
                                         entry represents an inequality:
                                         '=','>','<','>=','<='
            num_vertices -- (Integer) The desired number of vertices.
                            (List) Format: [<String>,<Integer>] WHERE the first entry
                            represents an inequality: '=','>','<','>=','<='
            perfect -- (Boolean)
            planar -- (Boolean)
            radius -- (Integer) The desired radius.
                      (List) Format: [<String>,<Integer>] WHERE the first entry represents
                             an inequality: '=','>','<','>=','<='
            regular -- (Boolean)
            spectrum -- (String) The desired spectrum.  (Ordered highest to lowest,
                                 delimited by ', ' and rounded to 6 decimal places).
            vertex_connectivity -- (Integer) The desired vertex connectivity.
                                   (List) Format: [<String>,<Integer>] WHERE the first
                                          entry represents an inequality:
                                          '=','>','<','>=','<='
            vertex_transitive -- (Boolean)

        EXAMPLES:
        TODO
        The basics:
            sage.: graphs_query.display_tables(tables=['graph_data','misc'], num_vertices=5,\
            ...                             lovasz_number=3.0, girth=4, radius=2, diameter=3)
            sage.: graphs_query.display_tables(tables=['degrees','spectrum','aut_grp','misc'], \
            ...                             layout='spring', num_hamiltonian_cycles=2,\
            ...                             regular=True, perfect=False)
            sage.: graphs_query.display_tables(tables=['degrees'],layout='spring',\
            ...                                degree_sequence=433211)

        Using Inequalities:
            sage.: graphs_query.display_tables(tables=['spectrum'], layout='circular', \
            ...                             min_eigenvalue=['=',-1], eigenvalues_sd=['<=',1], \
            ...                             energy=['>',5])

        The query string:
            sage.: graphs_query.display_tables(tables=['graph_data'], layout='spring', \
            ...             query='SELECT graph_data.graph6 \
            ...             FROM graph_data WHERE num_vertices<=4 \
            ...             and num_edges>3')
            sage.: graphs_query.display_tables(query='SELECT graph_data.graph6 FROM graph_data \
            ...             INNER JOIN degrees on graph_data.graph_id=degrees.graph_id \
            ...             WHERE num_vertices>6 and eulerian=1 and regular=0 and planar=1 \
            ...             and num_cycles<=2', tables=['graph_data','degrees'])
            sage.: graphs_query.display_tables(query="SELECT graph_data.graph6 \
            ...             FROM graph_data INNER JOIN misc on \
            ...             misc.graph_id=graph_data.graph_id WHERE \
            ...             misc.induced_subgraphs regexp '.*E~~w.*'", \
            ...             tables=['graph_data','misc','spectrum','degrees','aut_grp'])
        """
        from sage.plot.plot import plot

        if ( query is None):
            param = None
            query = __query_string__(graph6=graph6, num_vertices=num_vertices, num_edges=num_edges, num_cycles=num_cycles, num_hamiltonian_cycles=num_hamiltonian_cycles, eulerian=eulerian, planar=planar, perfect=perfect, lovasz_number=lovasz_number, complement_graph6=complement_graph6, aut_grp_size=aut_grp_size, num_orbits=num_orbits, num_fixed_points=num_fixed_points, vertex_transitive=vertex_transitive, edge_transitive=edge_transitive, degree_sequence=degree_sequence, min_degree=min_degree, max_degree=max_degree, average_degree=average_degree, degrees_sd=degrees_sd, regular=regular, vertex_connectivity=vertex_connectivity, edge_connectivity=edge_connectivity, num_components=num_components, girth=girth, radius=radius, diameter=diameter, clique_number=clique_number, independence_number=independence_number, num_cut_vertices=num_cut_vertices, min_vertex_cover_size=min_vertex_cover_size, num_spanning_trees=num_spanning_trees, induced_subgraphs=induced_subgraphs, spectrum=spectrum, min_eigenvalue=min_eigenvalue, max_eigenvalue=max_eigenvalue, eigenvalues_sd=eigenvalues_sd, energy=energy)
            query = re.sub('INNER JOIN .* WHERE','INNER JOIN aut_grp on aut_grp.graph_id=graph_data.graph_id INNER JOIN degrees on degrees.graph_id=graph_data.graph_id INNER JOIN misc on misc.graph_id=graph_data.graph_id INNER JOIN spectrum on spectrum.graph_id=graph_data.graph_id WHERE',query)
            query = re.sub('FROM graph_data WHERE','FROM graph_data INNER JOIN aut_grp on aut_grp.graph_id=graph_data.graph_id INNER JOIN degrees on degrees.graph_id=graph_data.graph_id INNER JOIN misc on misc.graph_id=graph_data.graph_id INNER JOIN spectrum on spectrum.graph_id=graph_data.graph_id WHERE',query)
            query = re.sub('SELECT graph_data.graph6','SELECT graph_data.graph6,graph_data.num_vertices,degrees.regular,aut_grp.aut_grp_size,graph_data.num_edges,degrees.min_degree,aut_grp.num_orbits,graph_data.num_cycles,degrees.max_degree,aut_grp.num_fixed_points,graph_data.num_hamiltonian_cycles,degrees.average_degree,aut_grp.vertex_transitive,graph_data.eulerian,degrees.degrees_sd,aut_grp.edge_transitive,graph_data.planar,degrees.degree_sequence,misc.vertex_connectivity,graph_data.perfect,spectrum.min_eigenvalue,misc.edge_connectivity,graph_data.lovasz_number,misc.girth,spectrum.max_eigenvalue,misc.num_cut_vertices,misc.independence_number,misc.radius,spectrum.eigenvalues_sd,misc.min_vertex_cover_size,misc.clique_number,misc.diameter,spectrum.energy,misc.num_spanning_trees,misc.num_components,graph_data.complement_graph6,spectrum.spectrum,misc.induced_subgraphs',query)
        else:
            # Deal only with the string:
            param = query.__param_tuple__
            query = query.__query_string__

        cur = (self.__connection__).cursor()
        if param is None:
            a = cur.execute(query)
            b = a.fetchall()
        else:
            tup = []
            # make it a tuple of strings:
            for i in range (len(param)):
                tup.append(str(param[i]))
            exe = cur.execute(query, tuple(tup))
            b = exe.fetchall()

        graph6list = []
        for i in range (len(b)):
            graph6 = str(b[i][0])
            g = graph.Graph('%s'%graph6)
            # The following line is time consuming and should not stay:
            graph6list.append(g.graph6_string())
            p = g.plot(layout=layout, vertex_size=30, vertex_labels=False, graph_border=False)
            p.save('%s.png'%i, figsize=[1,1])

        print "<html>"
        print '<table bgcolor=lightgrey cellpadding=0>'

        rows = 0
        for j in range(len(tables)):
            if ( tables[j] == 'misc' ): rows += 5
            elif ( tables[j] == 'graph_data' ): rows += 3
            elif ( tables[j] == 'aut_grp' or tables[j] == 'degrees' or tables[j] == 'spectrum' ): rows += 2

        from sage.misc.multireplace import multiple_replace
        to_bool = {'0':"False", '1':"True"}

        for i in range(len(b)):
            eul = multiple_replace(to_bool,'%s'%b[i][13])
            reg = multiple_replace(to_bool,'%s'%b[i][2])
            plan = multiple_replace(to_bool,'%s'%b[i][16])
            perf = multiple_replace(to_bool,'%s'%b[i][19])
            vtran = multiple_replace(to_bool,'%s'%b[i][12])
            etran = multiple_replace(to_bool,'%s'%b[i][15])

            print '<tr><td bgcolor=white align=center rowspan="%d"><img src="cell://%s.png"><br>%s</td>'%(rows,i,graph6list[i])
            top_row = True
            for j in range (len(tables)):
                if ( tables[j] == 'aut_grp' ):
                    if ( not top_row ): print '<tr>'
                    print '<td bgcolor=white align=left><font color=black> Aut Group Size: %s </font></td>'%(b[i][3])
                    print '<td bgcolor=white align=left><font color=black> Orbits: %s </font></td>'%(b[i][6])
                    print '<td bgcolor=white align=left><font color=black> Fixed Points: %s </font></td>'%(b[i][9])
                    print '</tr><tr>'
                    print '<td bgcolor=white align=left><font color=black> Vertex Transitive: %s </font></td>'%vtran
                    print '<td bgcolor=white align=left><font color=black> Edge Transitive: %s </font></td><td bgcolor=white></td></tr>'%etran
                    top_row = False
                elif ( tables[j] == 'degrees' ):
                    if ( not top_row ): print '<tr>'
                    print '<td bgcolor=white align=left><font color=black> Regular: %s </font></td>'%reg
                    print '<td bgcolor=white align=left><font color=black> Min Degree: %s </font></td>'%(b[i][5])
                    print '<td bgcolor=white align=left><font color=black> Average Degree: %s </font></td>'%(b[i][11])
                    print '</tr><tr>'
                    print '<td bgcolor=white align=left><font color=black> Degree Sequence: %s </font></td>'%(b[i][17])
                    print '<td bgcolor=white align=left><font color=black> Max Degree: %s </font></td>'%(b[i][8])
                    print '<td bgcolor=white align=left><font color=black> Degree SD: %s </font></td></tr>'%(b[i][14])
                    top_row = False
                elif ( tables[j] == 'graph_data' ):
                    if ( not top_row ): print '<tr>'
                    print '<td bgcolor=white align=left><font color=black> Vertices: %s </font></td>'%(b[i][1])
                    print '<td bgcolor=white align=left><font color=black> Hamiltonian Cycles: %s </font></td>'%(b[i][10])
                    print '<td bgcolor=white align=left><font color=black> Eulerian: %s </font></td>'%eul
                    print '</tr><tr>'
                    print '<td bgcolor=white align=left><font color=black> Edges: %s </font></td>'%(b[i][4])
                    print '<td bgcolor=white align=left><font color=black> Lovasz Number: %s </font></td>'%(b[i][22])
                    print '<td bgcolor=white align=left><font color=black> Planar: %s </font></td>'%plan
                    print '</tr><tr>'
                    print '<td bgcolor=white align=left><font color=black> Cycles: %s </font></td>'%(b[i][7])
                    print '<td bgcolor=white align=left><font color=black> Complement: %s </font></td>'%(b[i][35])
                    print '<td bgcolor=white align=left><font color=black> Perfect: %s </font></td></tr>'%perf
                    top_row = False
                elif ( tables[j] == 'misc' ):
                    if ( top_row == False ): print '<tr>'
                    print '<td bgcolor=white align=left><font color=black> Girth: %s </font></td>'%(b[i][23])
                    print '<td bgcolor=white align=left><font color=black> Clique Number: %s </font></td>'%(b[i][30])
                    print '<td bgcolor=white align=left><font color=black> Cut Vertices: %s </font></td>'%(b[i][25])
                    print '</tr><tr>'
                    print '<td bgcolor=white align=left><font color=black> Radius: %s </font></td>'%(b[i][27])
                    print '<td bgcolor=white align=left><font color=black> Independence Number: %s </font></td>'%(b[i][26])
                    print '<td bgcolor=white align=left><font color=black> Min Vertex Cover Size: %s </font></td>'%(b[i][29])
                    print '</tr><tr>'
                    print '<td bgcolor=white align=left><font color=black> Diameter: %s </font></td>'%(b[i][31])
                    print '<td bgcolor=white align=left><font color=black> Vertex Connectivity: %s </font></td>'%(b[i][18])
                    print '<td bgcolor=white align=left><font color=black> Spanning Trees: %s </font></td>'%(b[i][33])
                    print '</tr><tr>'
                    print '<td bgcolor=white align=left><font color=black> Components: %s </font></td>'%(b[i][34])
                    print '<td bgcolor=white align=left><font color=black> Edge Connectivity: %s </font></td><td bgcolor=white></td>'%(b[i][21])
                    print '</tr><tr><td bgcolor=white align=left colspan="3"><font color=black> Induced Subgraphs: %s </font></td></tr>'%(b[i][37])
                    top_row = False
                if ( tables[j] == 'spectrum' ):
                    if ( not top_row ): print '<tr>'
                    print '<td bgcolor=white align=left><font color=black> Min Eigenvalue: %s </font></td>'%(b[i][20])
                    print '<td bgcolor=white align=left><font color=black> Max Eigenvalue: %s </font></td>'%(b[i][24])
                    print '<td bgcolor=white align=left><font color=black> Eigenvalues SD: %s </font></td>'%(b[i][28])
                    print '</tr><tr>'
                    print '<td bgcolor=white align=left><font color=black> Energy: %s </font></td>'%(b[i][32])
                    print '<td bgcolor=white align=left colspan="2"><font color=black> Spectrum: %s </font></td></tr>'%(b[i][36])
                    top_row = False

            if ( i != len(b)-1 ): print '<tr><td bgcolor=lightblue colspan="4" height="3"></td></tr>'

        print "</table></html>"
        return

    def display_properties(self, properties=None, layout='circular', query=None, graph6=None, \
                           num_vertices=None, num_edges=None, num_cycles=None, num_hamiltonian_cycles=None, \
                           eulerian=None, planar=None, perfect=None, lovasz_number=None, \
                           complement_graph6=None, aut_grp_size=None, num_orbits=None, \
                           num_fixed_points=None, vertex_transitive=None, edge_transitive=None, \
                           degree_sequence=None, min_degree=None, max_degree=None, \
                           average_degree=None, degrees_sd=None, regular=None, \
                           vertex_connectivity=None, edge_connectivity=None, \
                           num_components=None, girth=None, radius=None, diameter=None, \
                           clique_number=None, independence_number=None, num_cut_vertices=None, \
                           min_vertex_cover_size=None, num_spanning_trees=None, \
                           induced_subgraphs=None, spectrum=None, min_eigenvalue=None, \
                           max_eigenvalue=None, eigenvalues_sd=None, energy=None):
        r"""
        Displays the results of a query in a table, including all specified
        properties and an image for each graph.

        INPUT:
            query -- (GenericSQLQuery) A sqlite query for graphs.db (See examples below).
            properties -- (List) A list of strings that are the exact name (as
                                 the following parameters) of the properties to
                                 display with the results.
            layout -- (String) The layout option for the graph image.  Options include:
                               'circular' -- plots the graph with vertices evenly
                                             distributed on a circle
                               'spring' -- uses the traditional spring layout
            aut_grp_size -- (Integer) The desired size of the automorphism group.
                            (List) Format: [<String>,<Integer>] WHERE the first
                                           entry represents an inequality:
                                           '=','>','<','>=','<='
            average_degree -- (Real) The desired average degree.
                              (List) Format: [<String>,<Real>] WHERE the first
                                             entry represents an inequality:
                                             '=','>','<','>=','<='
            clique_number -- (Integer) The desired clique number.
                             (List) Format: [<String>,<Integer>] WHERE the first
                                            entry represents an inequality:
                                            '=','>','<','>=','<='
            complement_graph6 -- (String) A graph6 string isomorphic to the
                                          desired complement graph.
                                 (List) A list of graph6 strings.  Will search
                                        for graphs with complement isomorphic to
                                        any string in the list.
            degree_sequence -- (Integer) The desired sequence of degrees.
                                         (Ordered highest to lowest).
            degrees_sd -- (Real) The desired standard deviation of degrees.
                          (List) Format: [<String>,<Real>] WHERE the first
                                         entry represents an inequality:
                                         '=','>','<','>=','<='
            diameter -- (Real) The desired diameter.
                        (List) Format: [<String>,<Real>] WHERE the first
                                        entry represents an inequality:
                                        '=','>','<','>=','<='
            edge_connectivity -- (Integer) The desired edge connectivity.
                                 (List) Format: [<String>,<Integer>] WHERE the first
                                        entry represents an inequality:
                                        '=','>','<','>=','<='
            edge_transitive -- (Boolean)
            eigenvalues_sd -- (Real) The desired standard deviation of eigenvalues.
                              (List) Format: [<String>,<Real>] WHERE the first
                                     entry represents an inequality:
                                     '=','>','<','>=','<='
            energy -- (Real) The desired energy.
                      (List) Format: [<String>,<Real>] WHERE the first entry
                                     represents an inequality: '=','>','<','>=','<='
            eulerian -- (Boolean)
            girth -- (Integer) The desired girth.
                     (List) Format: [<String>,<Integer>] WHERE the first entry
                                    represents an inequality: '=','>','<','>=','<='
            graph6 -- (String) A graph6 string isomorphic to the desired graph.
                      (List) A list of graph6 strings.  Will search for graphs
                             isomorphic to any string in the list.
            independence_number -- (Integer) The desired independence number.
                                   (List) Format: [<String>,<Integer>] WHERE the
                                   first entry represents an inequality:
                                   '=','>','<','>=','<='
            induced_subgraphs -- (String) graph6 string isomorphic to desired subgraph.
                                 (List) Format options:
                                        1. ['one_of',<String>,...,<String>]
                                           Will search for graphs containing a subgraph
                                           isomorphic to any of the graph6 strings in
                                           the list.
                                        2. ['all_of',<String>,...,<String>]
                                           Will search for graphs containing a subgraph
                                           isomorphic to each of the graph6 strings in
                                           the list.
            lovasz_number -- (Real) The desired lovasz number.
                             (List) Format: [<String>,<Real>] WHERE the first entry
                                    represents an inequality: '=','>','<','>=','<='
            max_degree -- (Integer) The desired maximum degree.
                          (List) Format: [<String>,<Integer>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
            max_eigenvalue -- (Real) The desired maximum eigenvalue.
                              (List) Format: [<String>,<Real>] WHERE the first entry
                                     represents an inequality: '=','>','<','>=','<='
            min_degree -- (Integer) The desired minimum degree.
                          (List) Format: [<String>,<Integer>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
            min_eigenvalue -- (Real) The desired minimum eigenvalue.
                              (List) Format: [<String>,<Real>] WHERE the first entry
                                     represents an inequality: '=','>','<','>=','<='
            min_vertex_cover_size -- (Integer) The desired minimum vertex cover size.
                                     (List) Format: [<String>,<Integer>] WHERE the
                                            first entry represents an inequality:
                                            '=','>','<','>=','<='
            num_components -- (Integer) The desired number of components.
                              (List) Format: [<String>,<Integer>] WHERE the first
                                     entry represents an inequality:
                                     '=','>','<','>=','<='
            num_cut_vertices -- (Integer) The desired number of cut vertices.
                                (List) Format: [<String>,<Integer>] WHERE the first
                                     entry represents an inequality:
                                     '=','>','<','>=','<='
            num_cycles -- (Integer) The desired number of cycles.
                          (List) Format: [<String>,<Integer>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
            num_edges -- (Integer) The desired number of edges.
                         (List) Format: [<String>,<Integer>] WHERE the first entry
                                represents an inequality: '=','>','<','>=','<='
            num_fixed_points -- (Integer) The desired number of fixed points.
                                (List) Format: [<String>,<Integer>] WHERE the first
                                       entry represents an inequality:
                                       '=','>','<','>=','<='
            num_hamiltonian_cycles -- (Integer) The desired number of hamiltonian cycles.
                                      (List) Format: [<String>,<Integer>] WHERE the first
                                             entry represents an inequality:
                                             '=','>','<','>=','<='
            num_orbits -- (Integer) The desired number of orbits.
                          (List) Format: [<String>,<Integer>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
            num_spanning_trees -- (Integer) The desired number of spanning trees.
                                  (List) Format: [<String>,<Integer>] WHERE the first
                                         entry represents an inequality:
                                         '=','>','<','>=','<='
            num_vertices -- (Integer) The desired number of vertices.
                            (List) Format: [<String>,<Integer>] WHERE the first entry
                            represents an inequality: '=','>','<','>=','<='
            perfect -- (Boolean)
            planar -- (Boolean)
            radius -- (Integer) The desired radius.
                      (List) Format: [<String>,<Integer>] WHERE the first entry represents
                             an inequality: '=','>','<','>=','<='
            regular -- (Boolean)
            spectrum -- (String) The desired spectrum.  (Ordered highest to lowest,
                                 delimited by ', ' and rounded to 6 decimal places).
            vertex_connectivity -- (Integer) The desired vertex connectivity.
                                   (List) Format: [<String>,<Integer>] WHERE the first
                                          entry represents an inequality:
                                          '=','>','<','>=','<='
            vertex_transitive -- (Boolean)

        EXAMPLES:
        TODO
        The basics:
            sage.: graphs_query.display_properties(properties=['num_vertices','lovasz_number',\
            ...                             'girth','radius','diameter'], num_vertices=5,\
            ...                             lovasz_number=3.0, girth=4, radius=2, diameter=3)
            sage.: graphs_query.display_properties(properties=['num_hamiltonian_cycles','regular',\
            ...                             'perfect','num_cycles','num_edges','spectrum'], \
            ...                             layout='spring', num_hamiltonian_cycles=2,\
            ...                             regular=True, perfect=False)
            sage.: graphs_query.display_properties(properties=['min_degree','max_degree',\
            ...                                'degrees_sd','average_degree','regular',\
            ...                                'induced_subgraphs'],layout='spring',\
            ...                                degree_sequence=433211)

        Using Inequalities:
            sage.: graphs_query.display_properties(properties=['energy','spectrum','eigenvalues_sd',\
            ...                             'complement_graph6'], layout='circular', \
            ...                             min_eigenvalue=['=',-1], eigenvalues_sd=['<=',1], \
            ...                             energy=['>',5])

        The query string:
            sage.: graphs_query.display_properties(properties=['eulerian','perfect','planar','regular',\
            ...             'edge_transitive','vertex_transitive','num_cycles','degree_sequence',\
            ...             'induced_subgraphs','num_vertices','max_degree'], layout='spring', \
            ...             query='SELECT graph_data.graph6 \
            ...             FROM graph_data WHERE num_vertices<=4 \
            ...             and num_edges>3')
            sage.: graphs_query.display_properties(query='SELECT graph_data.graph6 FROM graph_data \
            ...             INNER JOIN degrees on graph_data.graph_id=degrees.graph_id \
            ...             WHERE num_vertices>6 and eulerian=1 and regular=0 and planar=1 \
            ...             and num_cycles<=2', properties=['clique_number','independence_number'])
            sage.: graphs_query.display_properties(query="SELECT graph_data.graph6 \
            ...             FROM graph_data INNER JOIN misc on \
            ...             misc.graph_id=graph_data.graph_id WHERE \
            ...             misc.induced_subgraphs regexp '.*E~~w.*'", \
            ...             properties=['induced_subgraphs'])
        """
        from sage.plot.plot import plot

        if ( query is None):
            param = None
            query = __query_string__(graph6=graph6, num_vertices=num_vertices, num_edges=num_edges, num_cycles=num_cycles, num_hamiltonian_cycles=num_hamiltonian_cycles, eulerian=eulerian, planar=planar, perfect=perfect, lovasz_number=lovasz_number, complement_graph6=complement_graph6, aut_grp_size=aut_grp_size, num_orbits=num_orbits, num_fixed_points=num_fixed_points, vertex_transitive=vertex_transitive, edge_transitive=edge_transitive, degree_sequence=degree_sequence, min_degree=min_degree, max_degree=max_degree, average_degree=average_degree, degrees_sd=degrees_sd, regular=regular, vertex_connectivity=vertex_connectivity, edge_connectivity=edge_connectivity, num_components=num_components, girth=girth, radius=radius, diameter=diameter, clique_number=clique_number, independence_number=independence_number, num_cut_vertices=num_cut_vertices, min_vertex_cover_size=min_vertex_cover_size, num_spanning_trees=num_spanning_trees, induced_subgraphs=induced_subgraphs, spectrum=spectrum, min_eigenvalue=min_eigenvalue, max_eigenvalue=max_eigenvalue, eigenvalues_sd=eigenvalues_sd, energy=energy)
            query = re.sub('INNER JOIN .* WHERE','INNER JOIN aut_grp on aut_grp.graph_id=graph_data.graph_id INNER JOIN degrees on degrees.graph_id=graph_data.graph_id INNER JOIN misc on misc.graph_id=graph_data.graph_id INNER JOIN spectrum on spectrum.graph_id=graph_data.graph_id WHERE',query)
            query = re.sub('FROM graph_data WHERE','FROM graph_data INNER JOIN aut_grp on aut_grp.graph_id=graph_data.graph_id INNER JOIN degrees on degrees.graph_id=graph_data.graph_id INNER JOIN misc on misc.graph_id=graph_data.graph_id INNER JOIN spectrum on spectrum.graph_id=graph_data.graph_id WHERE',query)
            query = re.sub('SELECT graph_data.graph6','SELECT graph_data.graph6,graph_data.num_vertices,degrees.regular,aut_grp.aut_grp_size,graph_data.num_edges,degrees.min_degree,aut_grp.num_orbits,graph_data.num_cycles,degrees.max_degree,aut_grp.num_fixed_points,graph_data.num_hamiltonian_cycles,degrees.average_degree,aut_grp.vertex_transitive,graph_data.eulerian,degrees.degrees_sd,aut_grp.edge_transitive,graph_data.planar,degrees.degree_sequence,misc.vertex_connectivity,graph_data.perfect,spectrum.min_eigenvalue,misc.edge_connectivity,graph_data.lovasz_number,misc.girth,spectrum.max_eigenvalue,misc.num_cut_vertices,misc.independence_number,misc.radius,spectrum.eigenvalues_sd,misc.min_vertex_cover_size,misc.clique_number,misc.diameter,spectrum.energy,misc.num_spanning_trees,misc.num_components,graph_data.complement_graph6,spectrum.spectrum,misc.induced_subgraphs',query)
        else:
            # Deal only with the string:
            param = query.__param_tuple__
            query = query.__query_string__

        cur = (self.__connection__).cursor()
        if param is None:
            a = cur.execute(query)
            b = a.fetchall()
        else:
            tup = []
            # make it a tuple of strings:
            for i in range (len(param)):
                tup.append(str(param[i]))
            exe = cur.execute(query, tuple(tup))
            b = exe.fetchall()

        graph6list = []
        for i in range (len(b)):
            graph6 = str(b[i][0])
            g = graph.Graph('%s'%graph6)
            # The following line is time consuming and should not stay:
            graph6list.append(g.graph6_string())
            p = g.plot(layout=layout, vertex_size=30, vertex_labels=False, graph_border=False)
            p.save('%s.png'%i, figsize=[1,1])

        cells = len(properties)
        rows = 0
        max_cols = 0

        for i in range(len(properties)):
            if ( properties[i] == 'spectrum' or properties[i] == 'induced_subgraphs' ):
                properties.append(properties[i])
                properties[i] = 'skip'
                cells -= 1
                rows += 1
            elif (max_cols < 3):
                max_cols += 1

        if (max_cols == 0): max_cols = 1
        rows += (cells)/3
        if ( cells%3 != 0 ): rows += 1

        print "<html>"
        print '<table bgcolor=lightgrey cellpadding=0>'

        from sage.misc.multireplace import multiple_replace
        to_bool = {'0':"False", '1':"True"}

        for i in range(len(b)):
            eul = multiple_replace(to_bool,'%s'%b[i][13])
            reg = multiple_replace(to_bool,'%s'%b[i][2])
            plan = multiple_replace(to_bool,'%s'%b[i][16])
            perf = multiple_replace(to_bool,'%s'%b[i][19])
            vtran = multiple_replace(to_bool,'%s'%b[i][12])
            etran = multiple_replace(to_bool,'%s'%b[i][15])

            print '<tr><td bgcolor=white align=center rowspan="%d"><img src="cell://%s.png"><br>%s</td>'%(rows,i,graph6list[i])

            top_row = True
            index = 0
            j = 0
            while ( j < rows ):
                if ( not top_row ): print '<tr>'
                top_row = False
                if ( properties[index] == 'spectrum' ):
                    print '<td bgcolor=white align=left colspan="%d"><font color=black> Spectrum: %s </font></td></tr>'%(max_cols,b[i][36])
                    index += 1
                    j += 1
                elif ( properties[index] == 'induced_subgraphs' ):
                    print '<td bgcolor=white align=left colspan="%d"><font color=black> Induced Subgraphs: %s </font></td></tr>'%(max_cols,b[i][37])
                    index += 1
                    j += 1
                else:
                    for k in range (max_cols):
                        while ( properties[index] == 'skip' ): index += 1
                        if ( properties[index] == 'spectrum' or properties[index] == 'induced_subgraphs' ):
                            if ( k == 0 ):
                                j -= 1
                                top_row = True
                                break
                            else:
                                for m in range (max_cols-k):
                                    print '<td bgcolor=white></td>'
                                break
                        print '<td bgcolor=white align=left><font color=black> '
                        if ( properties[index] == 'num_vertices' ): print 'Vertices: %s '%b[i][1]
                        elif ( properties[index] == 'regular' ): print 'Regular: %s '%reg
                        elif ( properties[index] == 'aut_grp_size' ): print 'Aut Group Size: %s '%(b[i][3])
                        elif ( properties[index] == 'num_edges' ): print 'Edges: %s '%(b[i][4])
                        elif ( properties[index] == 'min_degree' ): print 'Min Degree: %s '%(b[i][5])
                        elif ( properties[index] == 'num_orbits' ): print 'Orbits: %s '%(b[i][6])
                        elif ( properties[index] == 'num_cycles' ): print 'Cycles: %s '%(b[i][7])
                        elif ( properties[index] == 'max_degree' ): print 'Max Degree: %s '%(b[i][8])
                        elif ( properties[index] == 'num_fixed_points' ): print 'Fixed Points: %s '%(b[i][9])
                        elif ( properties[index] == 'num_hamiltonian_cycles' ): print 'Hamiltonian Cycles: %s '%(b[i][10])
                        elif ( properties[index] == 'average_degree' ): print 'Average Degree: %s '%(b[i][11])
                        elif ( properties[index] == 'vertex_transitive' ): print 'Vertex Transitive: %s '%vtran
                        elif ( properties[index] == 'eulerian' ): print 'Eulerian: %s '%eul
                        elif ( properties[index] == 'degrees_sd' ): print 'Degree SD: %s '%(b[i][14])
                        elif ( properties[index] == 'edge_transitive' ): print 'Edge Transitive: %s '%etran
                        elif ( properties[index] == 'planar' ): print 'Planar: %s '%plan
                        elif ( properties[index] == 'degree_sequence' ): print 'Degree Sequence: %s '%(b[i][17])
                        elif ( properties[index] == 'vertex_connectivity' ): print 'Vertex Connectivity: %s '%(b[i][18])
                        elif ( properties[index] == 'perfect' ): print 'Perfect: %s '%perf
                        elif ( properties[index] == 'min_eigenvalue' ): print 'Min Eigenvalue: %s '%(b[i][20])
                        elif ( properties[index] == 'edge_connectivity' ): print 'Edge Connectivity: %s '%(b[i][21])
                        elif ( properties[index] == 'lovasz_number' ): print 'Lovasz Number: %s '%(b[i][22])
                        elif ( properties[index] == 'girth' ): print 'Girth: %s '%(b[i][23])
                        elif ( properties[index] == 'max_eigenvalue' ): print 'Max Eigenvalue: %s '%(b[i][24])
                        elif ( properties[index] == 'num_cut_vertices' ): print 'Cut Vertices: %s '%(b[i][25])
                        elif ( properties[index] == 'independence_number' ): print 'Independence Number: %s '%(b[i][26])
                        elif ( properties[index] == 'radius' ): print 'Radius: %s '%(b[i][27])
                        elif ( properties[index] == 'eigenvalues_sd' ): print 'Eigenvalues SD: %s '%(b[i][28])
                        elif ( properties[index] == 'min_vertex_cover_size' ): print 'Min Vertex Cover Size: %s '%(b[i][29])
                        elif ( properties[index] == 'clique_number' ): print 'Clique Number: %s '%(b[i][30])
                        elif ( properties[index] == 'diameter' ): print 'Diameter: %s '%(b[i][31])
                        elif ( properties[index] == 'energy' ): print 'Energy: %s '%(b[i][32])
                        elif ( properties[index] == 'num_spanning_trees' ): print 'Spanning Trees: %s '%(b[i][33])
                        elif ( properties[index] == 'num_components' ): print 'Components: %s '%(b[i][34])
                        elif ( properties[index] == 'complement_graph6' ): print 'Complement: %s '%(b[i][35])
                        print '</font></td>'
                        index += 1
                        if ( index >= len(properties) ):
                            if ( k != max_cols-1 ):
                                for m in range (max_cols-k):
                                    print '<td bgcolor=white></td>'
                            break
                    if ( not top_row ):
                        print '</tr>'
                        top_row = False
                    j += 1

            if ( i != len(b)-1 ): print '<tr><td bgcolor=lightblue colspan="%d" height="3"></td></tr>'%(max_cols+1)

        print "</table></html>"
        return

    def get_list(self, query=None, graph6=None, num_vertices=None, num_edges=None, \
                 num_cycles=None, num_hamiltonian_cycles=None, eulerian=None, planar=None, \
                 perfect=None, lovasz_number=None, complement_graph6=None, aut_grp_size=None, \
                 num_orbits=None, num_fixed_points=None, vertex_transitive=None, \
                 edge_transitive=None, degree_sequence=None, min_degree=None, \
                 max_degree=None, average_degree=None, degrees_sd=None, regular=None, \
                 vertex_connectivity=None, edge_connectivity=None, num_components=None, \
                 girth=None, radius=None, diameter=None, clique_number=None, \
                 independence_number=None, num_cut_vertices=None, min_vertex_cover_size=None, \
                 num_spanning_trees=None, induced_subgraphs=None, spectrum=None, \
                 min_eigenvalue=None, max_eigenvalue=None, eigenvalues_sd=None, energy=None):
        r"""
        Returns a list of SAGE graphs according to provided parameters.

        INPUT:
            query -- (GenericSQLQuery) A sqlite query for graphs.db (See examples below).
            aut_grp_size -- (Integer) The desired size of the automorphism group.
                            (List) Format: [<String>,<Integer>] WHERE the first
                                           entry represents an inequality:
                                           '=','>','<','>=','<='
            average_degree -- (Real) The desired average degree.
                              (List) Format: [<String>,<Real>] WHERE the first
                                             entry represents an inequality:
                                             '=','>','<','>=','<='
            clique_number -- (Integer) The desired clique number.
                             (List) Format: [<String>,<Integer>] WHERE the first
                                            entry represents an inequality:
                                            '=','>','<','>=','<='
            complement_graph6 -- (String) A graph6 string isomorphic to the
                                          desired complement graph.
                                 (List) A list of graph6 strings.  Will search
                                        for graphs with complement isomorphic to
                                        any string in the list.
            degree_sequence -- (Integer) The desired sequence of degrees.
                                         (Ordered highest to lowest).
            degrees_sd -- (Real) The desired standard deviation of degrees.
                          (List) Format: [<String>,<Real>] WHERE the first
                                         entry represents an inequality:
                                         '=','>','<','>=','<='
            diameter -- (Real) The desired diameter.
                        (List) Format: [<String>,<Real>] WHERE the first
                                        entry represents an inequality:
                                        '=','>','<','>=','<='
            edge_connectivity -- (Integer) The desired edge connectivity.
                                 (List) Format: [<String>,<Integer>] WHERE the first
                                        entry represents an inequality:
                                        '=','>','<','>=','<='
            edge_transitive -- (Boolean)
            eigenvalues_sd -- (Real) The desired standard deviation of eigenvalues.
                              (List) Format: [<String>,<Real>] WHERE the first
                                     entry represents an inequality:
                                     '=','>','<','>=','<='
            energy -- (Real) The desired energy.
                      (List) Format: [<String>,<Real>] WHERE the first entry
                                     represents an inequality: '=','>','<','>=','<='
            eulerian -- (Boolean)
            girth -- (Integer) The desired girth.
                     (List) Format: [<String>,<Integer>] WHERE the first entry
                                    represents an inequality: '=','>','<','>=','<='
            graph6 -- (String) A graph6 string isomorphic to the desired graph.
                      (List) A list of graph6 strings.  Will search for graphs
                             isomorphic to any string in the list.
            independence_number -- (Integer) The desired independence number.
                                   (List) Format: [<String>,<Integer>] WHERE the
                                   first entry represents an inequality:
                                   '=','>','<','>=','<='
            induced_subgraphs -- (String) graph6 string isomorphic to desired subgraph.
                                 (List) Format options:
                                        1. ['one_of',<String>,...,<String>]
                                           Will search for graphs containing a subgraph
                                           isomorphic to any of the graph6 strings in
                                           the list.
                                        2. ['all_of',<String>,...,<String>]
                                           Will search for graphs containing a subgraph
                                           isomorphic to each of the graph6 strings in
                                           the list.
            lovasz_number -- (Real) The desired lovasz number.
                             (List) Format: [<String>,<Real>] WHERE the first entry
                                    represents an inequality: '=','>','<','>=','<='
            max_degree -- (Integer) The desired maximum degree.
                          (List) Format: [<String>,<Integer>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
            max_eigenvalue -- (Real) The desired maximum eigenvalue.
                              (List) Format: [<String>,<Real>] WHERE the first entry
                                     represents an inequality: '=','>','<','>=','<='
            min_degree -- (Integer) The desired minimum degree.
                          (List) Format: [<String>,<Integer>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
            min_eigenvalue -- (Real) The desired minimum eigenvalue.
                              (List) Format: [<String>,<Real>] WHERE the first entry
                                     represents an inequality: '=','>','<','>=','<='
            min_vertex_cover_size -- (Integer) The desired minimum vertex cover size.
                                     (List) Format: [<String>,<Integer>] WHERE the
                                            first entry represents an inequality:
                                            '=','>','<','>=','<='
            num_components -- (Integer) The desired number of components.
                              (List) Format: [<String>,<Integer>] WHERE the first
                                     entry represents an inequality:
                                     '=','>','<','>=','<='
            num_cut_vertices -- (Integer) The desired number of cut vertices.
                                (List) Format: [<String>,<Integer>] WHERE the first
                                     entry represents an inequality:
                                     '=','>','<','>=','<='
            num_cycles -- (Integer) The desired number of cycles.
                          (List) Format: [<String>,<Integer>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
            num_edges -- (Integer) The desired number of edges.
                         (List) Format: [<String>,<Integer>] WHERE the first entry
                                represents an inequality: '=','>','<','>=','<='
            num_fixed_points -- (Integer) The desired number of fixed points.
                                (List) Format: [<String>,<Integer>] WHERE the first
                                       entry represents an inequality:
                                       '=','>','<','>=','<='
            num_hamiltonian_cycles -- (Integer) The desired number of hamiltonian cycles.
                                      (List) Format: [<String>,<Integer>] WHERE the first
                                             entry represents an inequality:
                                             '=','>','<','>=','<='
            num_orbits -- (Integer) The desired number of orbits.
                          (List) Format: [<String>,<Integer>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
            num_spanning_trees -- (Integer) The desired number of spanning trees.
                                  (List) Format: [<String>,<Integer>] WHERE the first
                                         entry represents an inequality:
                                         '=','>','<','>=','<='
            num_vertices -- (Integer) The desired number of vertices.
                            (List) Format: [<String>,<Integer>] WHERE the first entry
                            represents an inequality: '=','>','<','>=','<='
            perfect -- (Boolean)
            planar -- (Boolean)
            radius -- (Integer) The desired radius.
                      (List) Format: [<String>,<Integer>] WHERE the first entry represents
                             an inequality: '=','>','<','>=','<='
            regular -- (Boolean)
            spectrum -- (String) The desired spectrum.  (Ordered highest to lowest,
                                 delimited by ', ' and rounded to 6 decimal places).
            vertex_connectivity -- (Integer) The desired vertex connectivity.
                                   (List) Format: [<String>,<Integer>] WHERE the first
                                          entry represents an inequality:
                                          '=','>','<','>=','<='
            vertex_transitive -- (Boolean)

        EXAMPLES:
        TODO
            sage: g = graphs_query.get_list(num_vertices=5,lovasz_number=3.0,\
            ...                         girth=4,radius=2,diameter=3)
            ...
            sage: len(g)
            1
            sage.: g[0].show(layout='circular',figsize=[2,2],vertex_size=0,graph_border=True)
            sage: g = graphs_query.get_list(degree_sequence=433211)
            sage: graphs_list.to_graph6(g)
            'E@NW\nEAMw\n'
            sage: g = graphs_query.get_list(min_eigenvalue=['=',-1], \
            ...              eigenvalues_sd=['<=',1], energy=['>',5])
            ...
            sage: len(g) == graphs_query.number_of(min_eigenvalue=['=',-1], \
            ...              eigenvalues_sd=['<=',1], energy=['>',5])
            True
            sage: g = graphs_query.get_list(num_hamiltonian_cycles=2,regular=True,perfect=False)
            sage.: graphs_list.show_graphs(g)
            sage: g = graphs_query.get_list(query='SELECT graph_data.graph6 \
            ...             FROM graph_data WHERE num_vertices<=4 \
            ...             and num_edges>3')
            ...
            sage: graphs_list.to_graph6(g)
            'CN\nC]\nC^\nC~\n'
            sage: g = graphs_query.get_list(query='SELECT graph_data.graph6 FROM graph_data \
            ...             INNER JOIN degrees on graph_data.graph_id=degrees.graph_id \
            ...             WHERE num_vertices>6 and eulerian=1 and regular=0 and planar=1 \
            ...             and num_cycles<=2')
            ...
            sage: for i in range(len(g)):
            ...    if (g[i].is_isomorphic(Graph('FJ?GW'))):
            ...        print g[i].graph6_string()
            F@LAG
            sage: (Graph('FJ?GW')).is_isomorphic(Graph('F@LAG'))
            True
            sage.: g = graphs_query.get_list(query="SELECT graph_data.graph6 \
            ...             FROM graph_data INNER JOIN misc on \
            ...             misc.graph_id=graph_data.graph_id WHERE \
            ...             misc.induced_subgraphs regexp '.*E~~w.*'")
            ...
            sage.: graphs_list.show_graphs(g)
            sage.: graphs_list.to_graph6(g)
            'E~~w\nFJ\\zw\nFJ\\~w\nFJ^~w\nFJ~~w\nFN~~w\nF^~~w\nF~~~w\n'
        """
        if ( query is None ):
            param = None
            query = __query_string__(graph6=graph6, num_vertices=num_vertices, num_edges=num_edges, num_cycles=num_cycles, num_hamiltonian_cycles=num_hamiltonian_cycles, eulerian=eulerian, planar=planar, perfect=perfect, lovasz_number=lovasz_number, complement_graph6=complement_graph6, aut_grp_size=aut_grp_size, num_orbits=num_orbits, num_fixed_points=num_fixed_points, vertex_transitive=vertex_transitive, edge_transitive=edge_transitive, degree_sequence=degree_sequence, min_degree=min_degree, max_degree=max_degree, average_degree=average_degree, degrees_sd=degrees_sd, regular=regular, vertex_connectivity=vertex_connectivity, edge_connectivity=edge_connectivity, num_components=num_components, girth=girth, radius=radius, diameter=diameter, clique_number=clique_number, independence_number=independence_number, num_cut_vertices=num_cut_vertices, min_vertex_cover_size=min_vertex_cover_size, num_spanning_trees=num_spanning_trees, induced_subgraphs=induced_subgraphs, spectrum=spectrum, min_eigenvalue=min_eigenvalue, max_eigenvalue=max_eigenvalue, eigenvalues_sd=eigenvalues_sd, energy=energy)
        else:
            param = query.__param_tuple__
            query = query.__query_string__

        cur = (self.__connection__).cursor()
        if param is None:
            a = cur.execute(query)
            b = a.fetchall()
        else:
            tup = []
            # make it a tuple of strings:
            for i in range (len(param)):
                tup.append(str(param[i]))
            exe = cur.execute(query, tuple(tup))
            b = exe.fetchall()

        glist = []
        for i in range (len(b)):
            glist.append(graph.Graph(str(b[i][0])))
        return glist

    def number_of(self, query=None, graph6=None, num_vertices=None, num_edges=None, \
                  num_cycles=None, num_hamiltonian_cycles=None, eulerian=None, planar=None, \
                  perfect=None, lovasz_number=None, complement_graph6=None, aut_grp_size=None, \
                  num_orbits=None, num_fixed_points=None, vertex_transitive=None, \
                  edge_transitive=None, degree_sequence=None, min_degree=None, max_degree=None, \
                  average_degree=None, degrees_sd=None, regular=None, vertex_connectivity=None, \
                  edge_connectivity=None, num_components=None, girth=None, radius=None, \
                  diameter=None, clique_number=None, independence_number=None, \
                  num_cut_vertices=None, min_vertex_cover_size=None, num_spanning_trees=None, \
                  induced_subgraphs=None, spectrum=None, min_eigenvalue=None, max_eigenvalue=None, \
                  eigenvalues_sd=None, energy=None):
        r"""
        Returns the integer that represents the number of unlabeled graphs with 7 or
        fewer vertices that meet the provided search parameters.

        INPUT:
            query -- (GenericSQLQuery) A sqlite query for graphs.db (See examples below).
            aut_grp_size -- (Integer) The desired size of the automorphism group.
                            (List) Format: [<String>,<Integer>] WHERE the first
                                           entry represents an inequality:
                                           '=','>','<','>=','<='
            average_degree -- (Real) The desired average degree.
                              (List) Format: [<String>,<Real>] WHERE the first
                                             entry represents an inequality:
                                             '=','>','<','>=','<='
            clique_number -- (Integer) The desired clique number.
                             (List) Format: [<String>,<Integer>] WHERE the first
                                            entry represents an inequality:
                                            '=','>','<','>=','<='
            complement_graph6 -- (String) A graph6 string isomorphic to the
                                          desired complement graph.
                                 (List) A list of graph6 strings.  Will search
                                        for graphs with complement isomorphic to
                                        any string in the list.
            degree_sequence -- (Integer) The desired sequence of degrees.
                                         (Ordered highest to lowest).
            degrees_sd -- (Real) The desired standard deviation of degrees.
                          (List) Format: [<String>,<Real>] WHERE the first
                                         entry represents an inequality:
                                         '=','>','<','>=','<='
            diameter -- (Real) The desired diameter.
                        (List) Format: [<String>,<Real>] WHERE the first
                                        entry represents an inequality:
                                        '=','>','<','>=','<='
            edge_connectivity -- (Integer) The desired edge connectivity.
                                 (List) Format: [<String>,<Integer>] WHERE the first
                                        entry represents an inequality:
                                        '=','>','<','>=','<='
            edge_transitive -- (Boolean)
            eigenvalues_sd -- (Real) The desired standard deviation of eigenvalues.
                              (List) Format: [<String>,<Real>] WHERE the first
                                     entry represents an inequality:
                                     '=','>','<','>=','<='
            energy -- (Real) The desired energy.
                      (List) Format: [<String>,<Real>] WHERE the first entry
                                     represents an inequality: '=','>','<','>=','<='
            eulerian -- (Boolean)
            girth -- (Integer) The desired girth.
                     (List) Format: [<String>,<Integer>] WHERE the first entry
                                    represents an inequality: '=','>','<','>=','<='
            graph6 -- (String) A graph6 string isomorphic to the desired graph.
                      (List) A list of graph6 strings.  Will search for graphs
                             isomorphic to any string in the list.
            independence_number -- (Integer) The desired independence number.
                                   (List) Format: [<String>,<Integer>] WHERE the
                                   first entry represents an inequality:
                                   '=','>','<','>=','<='
            induced_subgraphs -- (String) graph6 string isomorphic to desired subgraph.
                                 (List) Format options:
                                        1. ['one_of',<String>,...,<String>]
                                           Will search for graphs containing a subgraph
                                           isomorphic to any of the graph6 strings in
                                           the list.
                                        2. ['all_of',<String>,...,<String>]
                                           Will search for graphs containing a subgraph
                                           isomorphic to each of the graph6 strings in
                                           the list.
            lovasz_number -- (Real) The desired lovasz number.
                             (List) Format: [<String>,<Real>] WHERE the first entry
                                    represents an inequality: '=','>','<','>=','<='
            max_degree -- (Integer) The desired maximum degree.
                          (List) Format: [<String>,<Integer>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
            max_eigenvalue -- (Real) The desired maximum eigenvalue.
                              (List) Format: [<String>,<Real>] WHERE the first entry
                                     represents an inequality: '=','>','<','>=','<='
            min_degree -- (Integer) The desired minimum degree.
                          (List) Format: [<String>,<Integer>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
            min_eigenvalue -- (Real) The desired minimum eigenvalue.
                              (List) Format: [<String>,<Real>] WHERE the first entry
                                     represents an inequality: '=','>','<','>=','<='
            min_vertex_cover_size -- (Integer) The desired minimum vertex cover size.
                                     (List) Format: [<String>,<Integer>] WHERE the
                                            first entry represents an inequality:
                                            '=','>','<','>=','<='
            num_components -- (Integer) The desired number of components.
                              (List) Format: [<String>,<Integer>] WHERE the first
                                     entry represents an inequality:
                                     '=','>','<','>=','<='
            num_cut_vertices -- (Integer) The desired number of cut vertices.
                                (List) Format: [<String>,<Integer>] WHERE the first
                                     entry represents an inequality:
                                     '=','>','<','>=','<='
            num_cycles -- (Integer) The desired number of cycles.
                          (List) Format: [<String>,<Integer>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
            num_edges -- (Integer) The desired number of edges.
                         (List) Format: [<String>,<Integer>] WHERE the first entry
                                represents an inequality: '=','>','<','>=','<='
            num_fixed_points -- (Integer) The desired number of fixed points.
                                (List) Format: [<String>,<Integer>] WHERE the first
                                       entry represents an inequality:
                                       '=','>','<','>=','<='
            num_hamiltonian_cycles -- (Integer) The desired number of hamiltonian cycles.
                                      (List) Format: [<String>,<Integer>] WHERE the first
                                             entry represents an inequality:
                                             '=','>','<','>=','<='
            num_orbits -- (Integer) The desired number of orbits.
                          (List) Format: [<String>,<Integer>] WHERE the first entry
                                 represents an inequality: '=','>','<','>=','<='
            num_spanning_trees -- (Integer) The desired number of spanning trees.
                                  (List) Format: [<String>,<Integer>] WHERE the first
                                         entry represents an inequality:
                                         '=','>','<','>=','<='
            num_vertices -- (Integer) The desired number of vertices.
                            (List) Format: [<String>,<Integer>] WHERE the first entry
                            represents an inequality: '=','>','<','>=','<='
            perfect -- (Boolean)
            planar -- (Boolean)
            radius -- (Integer) The desired radius.
                      (List) Format: [<String>,<Integer>] WHERE the first entry represents
                             an inequality: '=','>','<','>=','<='
            regular -- (Boolean)
            spectrum -- (String) The desired spectrum.  (Ordered highest to lowest,
                                 delimited by ', ' and rounded to 6 decimal places).
            vertex_connectivity -- (Integer) The desired vertex connectivity.
                                   (List) Format: [<String>,<Integer>] WHERE the first
                                          entry represents an inequality:
                                          '=','>','<','>=','<='
            vertex_transitive -- (Boolean)

        EXAMPLES:
        TODO
            sage: graphs_query.number_of()
            1252
            sage: g = graphs_query.get_list(num_vertices=5,lovasz_number=3.0,\
            ...                         girth=4,radius=2,diameter=3)
            ...
            sage: h = graphs_query.number_of(num_vertices=5,lovasz_number=3.0,\
            ...                         girth=4,radius=2,diameter=3)
            ...
            sage: h == len(g)
            True
            sage: for i in range(8)[1:]:
            ...    g = graphs_query.number_of(num_vertices=i)
            ...    h = graphs_query.get_list(num_vertices=i)
            ...    print g == len(h)
            True
            True
            True
            True
            True
            True
            True
            sage: graphs_query.number_of(degree_sequence=433211)
            2
            sage: graphs_query.number_of(min_eigenvalue=['=',-1], \
            ...              eigenvalues_sd=['<=',1], energy=['>',5])
            2
            sage: graphs_query.number_of(num_hamiltonian_cycles=2,regular=True,perfect=False)
            2
            sage: graphs_query.number_of(query='SELECT graph_data.graph6 \
            ...             FROM graph_data WHERE num_vertices<=4 \
            ...             and num_edges>3')
            4
            sage: graphs_query.number_of(query='SELECT graph_data.graph6 FROM graph_data \
            ...             INNER JOIN degrees on graph_data.graph_id=degrees.graph_id \
            ...             WHERE num_vertices>6 and eulerian=1 and regular=0 and planar=1 \
            ...             and num_cycles<=2')
            9
            sage: graphs_query.number_of(num_vertices=7,num_cycles=['>',2])
            913
            sage: a = graphs_query.number_of(num_hamiltonian_cycles=['>=',2])
            sage: b = graphs_query.number_of(num_hamiltonian_cycles=['<',2])
            sage: c = graphs_query.number_of(num_hamiltonian_cycles=['<=',40])
            sage: d = graphs_query.number_of(num_hamiltonian_cycles=['>',40])
            sage: a-d == c-b
            True
            sage: q = 'SELECT graph_data.graph6 FROM graph_data WHERE \
            ...                      num_hamiltonian_cycles=2 or '
            ...
            sage: for i in range(40)[3:]:
            ...    q += 'num_hamiltonian_cycles=%d or '%i
            ...
            sage: q += 'num_hamiltonian_cycles=40'
            sage: long = graphs_query.number_of(query=q)
            sage: short = graphs_query.number_of(query='SELECT graph_data.graph6 \
            ...                      FROM graph_data WHERE num_hamiltonian_cycles>=2 \
            ...                      and num_hamiltonian_cycles<=40')
            ...
            sage: long == short == a-d
            True
            sage: graphs_query.number_of(query="SELECT graph_data.graph6 \
            ...             FROM graph_data INNER JOIN misc on \
            ...             misc.graph_id=graph_data.graph_id WHERE \
            ...             misc.induced_subgraphs regexp '.*E~~w.*'")
            8
        """
        glist = self.get_list(query=query, graph6=graph6, num_vertices=num_vertices, num_edges=num_edges, num_cycles=num_cycles, num_hamiltonian_cycles=num_hamiltonian_cycles, eulerian=eulerian, planar=planar, perfect=perfect, lovasz_number=lovasz_number, complement_graph6=complement_graph6, aut_grp_size=aut_grp_size, num_orbits=num_orbits, num_fixed_points=num_fixed_points, vertex_transitive=vertex_transitive, edge_transitive=edge_transitive, degree_sequence=degree_sequence, min_degree=min_degree, max_degree=max_degree, average_degree=average_degree, degrees_sd=degrees_sd, regular=regular, vertex_connectivity=vertex_connectivity, edge_connectivity=edge_connectivity, num_components=num_components, girth=girth, radius=radius, diameter=diameter, clique_number=clique_number, independence_number=independence_number, num_cut_vertices=num_cut_vertices, min_vertex_cover_size=min_vertex_cover_size, num_spanning_trees=num_spanning_trees, induced_subgraphs=induced_subgraphs, spectrum=spectrum, min_eigenvalue=min_eigenvalue, max_eigenvalue=max_eigenvalue, eigenvalues_sd=eigenvalues_sd, energy=energy)
        return len(glist)

def __query_string__(graph6=None, num_vertices=None, num_edges=None, num_cycles=None, num_hamiltonian_cycles=None, eulerian=None, planar=None, perfect=None, lovasz_number=None, complement_graph6=None, aut_grp_size=None, num_orbits=None, num_fixed_points=None, vertex_transitive=None, edge_transitive=None, degree_sequence=None, min_degree=None, max_degree=None, average_degree=None, degrees_sd=None, regular=None, vertex_connectivity=None, edge_connectivity=None, num_components=None, girth=None, radius=None, diameter=None, clique_number=None, independence_number=None, num_cut_vertices=None, min_vertex_cover_size=None, num_spanning_trees=None, induced_subgraphs=None, spectrum=None, min_eigenvalue=None, max_eigenvalue=None, eigenvalues_sd=None, energy=None):
    """
    Creates the query string for sqlite database.  Will query graph_data automatically,
    and any other tables as necessary.

    Applies parameters specified in get_list, number_of and display functions.
    """
    query = 'SELECT graph_data.graph6 FROM graph_data WHERE graph_data.graph_id=graph_data.graph_id and'
    if ( aut_grp_size != None or num_orbits != None or num_fixed_points != None or vertex_transitive != None or edge_transitive != None ):
        query = __aut_grp_string__(query=query, aut_grp_size=aut_grp_size, num_orbits=num_orbits, num_fixed_points=num_fixed_points, vertex_transitive=vertex_transitive, edge_transitive=edge_transitive)
    if ( degree_sequence != None or min_degree != None or max_degree != None or average_degree != None or degrees_sd != None or regular != None ):
        query = __degrees_string__(query=query, degree_sequence=degree_sequence, min_degree=min_degree, max_degree=max_degree, average_degree=average_degree, degrees_sd=degrees_sd, regular=regular)
    if ( vertex_connectivity != None or edge_connectivity != None or num_components != None or girth != None or radius != None or diameter != None or clique_number != None or independence_number != None or num_cut_vertices != None or min_vertex_cover_size != None or num_spanning_trees != None or induced_subgraphs != None ):
        query = __misc_string__(query=query, vertex_connectivity=vertex_connectivity, edge_connectivity=edge_connectivity, num_components=num_components, girth=girth, radius=radius, diameter=diameter, clique_number=clique_number, independence_number=independence_number, num_cut_vertices=num_cut_vertices, min_vertex_cover_size=min_vertex_cover_size, num_spanning_trees=num_spanning_trees, subgraph=induced_subgraphs)
    if ( spectrum != None or min_eigenvalue != None or max_eigenvalue != None or eigenvalues_sd != None or energy != None ):
        query = __spectrum_string__(query=query, spectrum=spectrum, min_eigenvalue=min_eigenvalue, max_eigenvalue=max_eigenvalue, eigenvalues_sd=eigenvalues_sd, energy=energy)
    if ( graph6 != None ):
        if (str(type(graph6)) == "<type 'list'>"):
            # only one_of
            for i in range (len(graph6)):
                graph6[i] = ((graph.Graph(graph6[i])).canonical_label()).graph6_string()
            query += ' ('
            for i in range (len(graph6)-1):
                query += "graph_data.graph6='%s' or "%graph6[i]
            query += "graph_data.graph6='%s') and"%graph6[len(graph6)-1]
        else:
            graph6 = ((graph.Graph(graph6)).canonical_label()).graph6_string()
            query += " graph_data.graph6='%s' and"%graph6
    if ( num_vertices != None ):
        if (str(type(num_vertices)) == "<type 'list'>"):
            query += ' graph_data.num_vertices%s%d and'%(num_vertices[0], num_vertices[1])
        else:
            query += ' graph_data.num_vertices=%d and'%num_vertices
    if ( num_edges != None ):
        if (str(type(num_edges)) == "<type 'list'>"):
            query += ' graph_data.num_edges%s%d and'%(num_edges[0], num_edges[1])
        else:
            query += ' graph_data.num_edges=%d and'%num_edges
    if ( num_cycles != None ):
        if (str(type(num_cycles)) == "<type 'list'>"):
            query += ' graph_data.num_cycles%s%d and'%(num_cycles[0], num_cycles[1])
        else:
            query += ' graph_data.num_cycles=%d and'%num_cycles
    if ( num_hamiltonian_cycles != None ):
        if (str(type(num_hamiltonian_cycles)) == "<type 'list'>"):
            query += ' graph_data.num_hamiltonian_cycles%s%d and'%(num_hamiltonian_cycles[0], num_hamiltonian_cycles[1])
        else:
            query += ' graph_data.num_hamiltonian_cycles=%d and'%num_hamiltonian_cycles
    if ( eulerian != None ):
        query += ' graph_data.eulerian=%d and'%eulerian
    if ( planar != None ):
        query += ' graph_data.planar=%d and'%planar
    if ( perfect != None ):
        query += ' graph_data.perfect=%d and'%perfect
    if ( lovasz_number != None ):
        if (str(type(lovasz_number)) == "<type 'list'>"):
            query += ' graph_data.lovasz_number%s%d and'%(lovasz_number[0], lovasz_number[1])
        else:
            query += ' graph_data.lovasz_number=%d and'%lovasz_number
    if ( complement_graph6 != None ):
        if (str(type(complement_graph6)) == "<type 'list'>"):
            # only one_of
            for i in range (len(complement_graph6)):
                complement_graph6[i] = ((graph.Graph(complement_graph6[i])).canonical_label()).graph6_string()
            query += ' ('
            for i in range (len(complement_graph6)-1):
                query += "graph_data.complement_graph6='%s' or "%complement_graph6[i]
            query += "graph_data.complement_graph6='%s') and"%complement_graph6[len(complement_graph6)-1]
        else:
            complement_graph6 = ((graph.Graph(complement_graph6)).canonical_label()).graph6_string()
            query += " graph_data.complement_graph6='%s' and"%complement_graph6

    clean_query = re.sub(' and$','',query)
    return clean_query

def __aut_grp_string__(query=None, aut_grp_size=None, num_orbits=None, num_fixed_points=None, vertex_transitive=None, edge_transitive=None):
    """
    Appends necessary info to query string to search for graphs with
    properties specified in the aut_grp database table.

    Applies parameters specified in get_list, number_of and display functions.
    """
    r = query.split('WHERE',1)
    s = 'INNER JOIN aut_grp on aut_grp.graph_id=graph_data.graph_id WHERE'
    query = s.join(r)

    if ( aut_grp_size is not None ):
        if (str(type(aut_grp_size)) == "<type 'list'>"):
            query += ' aut_grp.aut_grp_size%s%d and'%(aut_grp_size[0], aut_grp_size[1])
        else:
            query += ' aut_grp.aut_grp_size=%d and'%aut_grp_size
    if ( num_orbits is not None ):
        if (str(type(num_orbits)) == "<type 'list'>"):
            query += ' aut_grp.num_orbits%s%d and'%(num_orbits[0], num_orbits[1])
        else:
            query += ' aut_grp.num_orbits=%d and'%num_orbits
    if ( num_fixed_points is not None ):
        if (str(type(num_fixed_points)) == "<type 'list'>"):
            query += ' aut_grp.num_fixed_points%s%d and'%(num_fixed_points[0], num_fixed_points[1])
        else:
            query += ' aut_grp.num_fixed_points=%d and'%num_fixed_points
    if ( vertex_transitive is not None ):
        query += ' aut_grp.vertex_transitive=%d and'%vertex_transitive
    if ( edge_transitive is not None ):
        query += ' aut_grp.edge_transitive=%d and'%edge_transitive

    return query

def __degrees_string__(query=None, degree_sequence=None, min_degree=None, max_degree=None, average_degree=None, degrees_sd=None, regular=None):
    """
    Appends necessary info to query string to search for graphs with
    properties specified in the degrees database table.

    Applies parameters specified in get_list, number_of and display functions.
    """
    r = query.split('WHERE',1)
    s = 'INNER JOIN degrees on degrees.graph_id=graph_data.graph_id WHERE'
    query = s.join(r)

    if ( degree_sequence is not None ):
        query += ' degrees.degree_sequence=%d and'%degree_sequence
    if ( min_degree is not None ):
        if (str(type(min_degree)) == "<type 'list'>"):
            query += ' degrees.min_degree%s%d and'%(min_degree[0], min_degree[1])
        else:
            query += ' degrees.min_degree=%d and'%min_degree
    if ( max_degree is not None ):
        if (str(type(max_degree)) == "<type 'list'>"):
            query += ' degrees.max_degree%s%d and'%(max_degree[0], max_degree[1])
        else:
            query += ' degrees.max_degree=%d and'%max_degree
    if ( average_degree is not None ):
        if (str(type(average_degree)) == "<type 'list'>"):
            query += ' degrees.average_degree%s%s and'%(average_degree[0], average_degree[1])
        else:
            query += ' degrees.average_degree=%s and'%average_degree
    if ( degrees_sd is not None ):
        if (str(type(degrees_sd)) == "<type 'list'>"):
            query += ' degrees.degrees_sd%s%s and'%(degrees_sd[0], degrees_sd[1])
        else:
            query += ' degrees.degrees_sd=%s and'%degrees_sd
    if ( regular is not None ):
        query += ' degrees.regular=%d and'%regular

    return query

def __misc_string__(query=None, vertex_connectivity=None, edge_connectivity=None, num_components=None, girth=None, radius=None, diameter=None, clique_number=None, independence_number=None, num_cut_vertices=None, min_vertex_cover_size=None, num_spanning_trees=None, subgraph=None):
    """
    Appends necessary info to query string to search for graphs with
    properties specified in the misc database table.

    Applies parameters specified in get_list, number_of and display functions.
    """
    r = query.split('WHERE',1)
    s = 'INNER JOIN misc on misc.graph_id=graph_data.graph_id WHERE'
    query = s.join(r)

    if (subgraph is not None):
        from sage.misc.multireplace import multiple_replace
        clean = {'[': '[[]', ']': '[]]', '?': '[?]', '{': '[{]', '}': '[}]', '^': '[^]', '|': '[|]'}
        if (str(type(subgraph)) == "<type 'list'>"):
            for i in range (len(subgraph))[1:]:
                subgraph[i] = ((graph.Graph(subgraph[i])).canonical_label()).graph6_string()
                subgraph[i] = multiple_replace(clean,subgraph[i])
                subgraph[i] = subgraph[i].replace('\\',"\\\\")
            if (subgraph[0] == 'one_of'):
                query += ' ('
                for i in range (len(subgraph)-1)[1:]:
                    query += "misc.induced_subgraphs regexp '.*%s.*' or "%subgraph[i]
                query += "misc.induced_subgraphs regexp '.*%s.*') and"%subgraph[len(subgraph)-1]
            elif (subgraph[0] == 'all_of'):
                for i in range (len(subgraph))[1:]:
                    query += " misc.induced_subgraphs regexp '.*%s.*' and"%subgraph[i]
        else:
            subgraph = ((graph.Graph(subgraph)).canonical_label()).graph6_string()
            subgraph = multiple_replace(clean,subgraph)
            subgraph = subgraph.replace('\\',"\\\\")
            query += " misc.induced_subgraphs regexp '.*%s.*' and"%subgraph

    if (vertex_connectivity is not None):
        if (str(type(vertex_connectivity)) == "<type 'list'>"):
            query += ' misc.vertex_connectivity%s%d and'%(vertex_connectivity[0], vertex_connectivity[1])
        else:
            query += ' misc.vertex_connectivity=%d and'%vertex_connectivity
    if (edge_connectivity is not None):
        if (str(type(edge_connectivity)) == "<type 'list'>"):
            query += ' misc.edge_connectivity%s%d and'%(edge_connectivity[0], edge_connectivity[1])
        else:
            query += ' misc.edge_connectivity=%d and'%edge_connectivity
    if (num_components is not None):
        if (str(type(num_components)) == "<type 'list'>"):
            query += ' misc.num_components%s%d and'%(num_components[0], num_components[1])
        else:
            query += ' misc.num_components=%d and'%num_components
    if (girth is not None):
        if (str(type(girth)) == "<type 'list'>"):
            query += ' misc.girth%s%d and'%(girth[0], girth[1])
        else:
            query += ' misc.girth=%d and'%girth
    if (radius is not None):
        if (str(type(radius)) == "<type 'list'>"):
            query += ' misc.radius%s%d and'%(radius[0], radius[1])
        else:
            query += ' misc.radius=%d and'%radius
    if (diameter is not None):
        if (str(type(diameter)) == "<type 'list'>"):
            query += ' misc.diameter%s%d and'%(diameter[0], diameter[1])
        else:
            query += ' misc.diameter=%d and'%diameter
    if (clique_number is not None):
        if (str(type(clique_number)) == "<type 'list'>"):
            query += ' misc.clique_number%s%d and'%(clique_number[0], clique_number[1])
        else:
            query += ' misc.clique_number=%d and'%clique_number
    if (independence_number is not None):
        if (str(type(independence_number)) == "<type 'list'>"):
            query += ' misc.independence_number%s%d and'%(independence_number[0], independence_number[1])
        else:
            query += ' misc.independence_number=%d and'%independence_number
    if (num_cut_vertices is not None):
        if (str(type(num_cut_vertices)) == "<type 'list'>"):
            query += ' misc.num_cut_vertices%s%d and'%(num_cut_vertices[0], num_cut_vertices[1])
        else:
            query += ' misc.num_cut_vertices=%d and'%num_cut_vertices
    if (min_vertex_cover_size is not None):
        if (str(type(min_vertex_cover_size)) == "<type 'list'>"):
            query += ' misc.min_vertex_cover_size%s%d and'%(min_vertex_cover_size[0], min_vertex_cover_size[1])
        else:
            query += ' misc.min_vertex_cover_size=%d and'%min_vertex_cover_size
    if (num_spanning_trees is not None):
        if (str(type(num_spanning_trees)) == "<type 'list'>"):
            query += ' misc.num_spanning_trees%s%d and'%(num_spanning_trees[0], num_spanning_trees[1])
        else:
            query += ' misc.num_spanning_trees=%d and'%num_spanning_trees

    return query

def __spectrum_string__(query=None, spectrum=None, min_eigenvalue=None, max_eigenvalue=None, eigenvalues_sd=None, energy=None):
    """
    Appends necessary info to query string to search for graphs with
    properties specified in the spectrum database table.

    Applies parameters specified in get_list, number_of and display functions.
    """
    r = query.split('WHERE',1)
    s = 'INNER JOIN spectrum on spectrum.graph_id=graph_data.graph_id WHERE'
    query = s.join(r)

    if ( spectrum is not None ):
        query += ' spectrum.spectrum=%s and'%spectrum
    if ( min_eigenvalue is not None ):
        if (str(type(min_eigenvalue)) == "<type 'list'>"):
            query += ' spectrum.min_eigenvalue%s%s and'%(min_eigenvalue[0], min_eigenvalue[1])
        else:
            query += ' spectrum.min_eigenvalue=%s and'%min_eigenvalue
    if ( max_eigenvalue is not None ):
        if (str(type(max_eigenvalue)) == "<type 'list'>"):
            query += ' spectrum.max_eigenvalue%s%s and'%(max_eigenvalue[0], max_eigenvalue[1])
        else:
            query += ' spectrum.max_eigenvalue=%s and'%max_eigenvalue
    if ( eigenvalues_sd is not None ):
        if (str(type(eigenvalues_sd)) == "<type 'list'>"):
            query += ' spectrum.eigenvalues_sd%s%s and'%(eigenvalues_sd[0], eigenvalues_sd[1])
        else:
            query += ' spectrum.eigenvalues_sd=%s and'%eigenvalues_sd
    if ( energy is not None ):
        if (str(type(energy)) == "<type 'list'>"):
            query += ' spectrum.energy%s%s and'%(energy[0], energy[1])
        else:
            query += ' spectrum.energy=%s and'%energy

    return query
