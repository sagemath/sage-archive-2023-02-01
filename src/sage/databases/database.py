"""
Relational (sqlite) Databases Module.

INFO:

    This module implements classes (SQLDatabase (mutable), GenericSQLDatabase
    (immutable), SQLQuery (pythonic implementation for the user with little or
    no knowledge of sqlite) and GenericSQLQuery (for the more advanced user))
    that wrap the basic functionality of sqlite.

    Databases are constructed via a triple indexed dictionary called a skeleton.
    A skeleton should be constructed to fit the following format:

    skeleton -- a triple-indexed dictionary
        outer key - table name
            inner key - column name
                inner inner key - one of the following:
                    primary_key - boolean, whether column has been set as primary key
                    index - boolean, whether column has been set as index
                    sql - one of 'TEXT', 'BOOLEAN', 'INTEGER', 'REAL', or other
                            user defined type

    An example skeleton of a database with one table, that table with one column:
        {'table1':{'col1':{'primary_key':False, 'index':True, 'sql':'REAL'}}}

    SQLDatabases can also be constructed via the add, drop and commit functions.
    The vacuum function is also useful for restoring hard disk space after a
    database has shrunk in size.

    A SQLQuery can be constructed by providing a query_dict, which is a dictionary
    with the following sample format:
        {'table_name': 'tblname', 'display_cols': ['col1', 'col2', 'col3'],
                                    'expression':[col, operator, value]}

    Finally a GenericSQLQuery allows the user to directly input the query string
    for a database, and also supports the '?' syntax by allowing an argument for
    a tuple of parameters to query.

    For full details, please see the tutorial.  Also, sage.graphs.graph_database.py
    is an example of implementing a database class in Sage using this interface.

AUTHORS:
    -- Emily A. Kirkman (2008-09-20): added functionality to generate plots and
                                    reformat output in show
    -- Emily A. Kirkman and Robert L. Miller (2007-06-17): initial version

"""
# FUTURE TODOs (Ignore for now):
#    - order by clause in query strings
#    - delete from query containing joins
#    - add data by column
#    - wrap sqlalchemy
#    - create query interface (with interact)
#    - allow kwds arguments to SQLQuery (like GraphQuery)

################################################################################
#           Copyright (C) 2007 Emily A. Kirkman
#                              Robert L. Miller
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
################################################################################
from sqlite3 import dbapi2 as sqlite # if anyone would like to explain why dbapi2...
import os
import re
from sage.misc.misc import tmp_filename
from sage.structure.sage_object import SageObject
import sage.server.support
from sage.plot.plot import plot
from sage.graphs.graph import Graph


def regexp(expr, item):
    """
    Function to define regular expressions in pysqlite.
    Returns 1 if parameter `item` matches the regular expression parameter `expr`.
    Returns 0 otherwise (i.e.: no match).

    REFERENCES:
        Gerhard Haring. [Online] Available: http://lists.initd.org/pipermail/pysqlite/2005-November/000253.html

    EXAMPLE:
        sage: from sage.databases.database import regexp
        sage: regexp('.s.*','cs')
        True
        sage: regexp('.s.*','ccs')
        False
        sage: regexp('.s.*','cscccc')
        True
    """
    r = re.compile(expr)
    return r.match(item) is not None

def verify_type(type):
    """
    Verify that the specified type is one of the allowed strings.

    EXAMPLES:
        sage: from sage.databases.database import verify_type
        sage: verify_type('INT')
        True
        sage: verify_type('int')
        True
        sage: verify_type('float')
        Traceback (most recent call last):
        ...
        TypeError: float is not a legal type.

    """
    types = ['INTEGER','INT','BOOLEAN','REAL','TEXT','BOOL','BLOB']
    if type.upper() not in types:
        raise TypeError('%s is not a legal type.'%type)
    return True

def verify_column(col_dict):
    """
    Verify that a column dict is in proper format*, and return a dict with
    default values filled in.

    * {'primary_key':False, 'index':False, 'sql':'REAL'}

    EXAMPLES:
        sage: from sage.databases.database import verify_column
        sage: col = {'sql':'BOOLEAN'}
        sage: verify_column(col)
        {'index': False, 'primary_key': False, 'sql': 'BOOLEAN'}
        sage: col = {'primary_key':True, 'sql':'INTEGER'}
        sage: verify_column(col)
        {'index': False, 'primary_key': True, 'sql': 'INTEGER'}
        sage: verify_column({})
        Traceback (most recent call last):
        ...
        ValueError: SQL type must be declared, e.g. {'sql':'REAL'}.

    """
    d = {}
    d['primary_key'] = col_dict.get('primary_key', False)
    d['index'] = col_dict.get('index', False)
    if not col_dict.has_key('sql'):
        raise ValueError("SQL type must be declared, e.g. {'sql':'REAL'}.")
    if verify_type(col_dict['sql']):
        d['sql'] = col_dict['sql']
    return d

def verify_operator(operator):
    """
    Checks that the provided operator is one of the allowed strings.
    Legal operators include the following strings:
        - '='
        - '<='
        - '>='
        - '<'
        - '>'
        - '<>'
        - 'like'
        - 'regexp'
        - 'is null'
        - 'is not null'

    EXAMPLES:
        sage: from sage.databases.database import verify_operator
        sage: verify_operator('<=')
        True
        sage: verify_operator('regexp')
        True
        sage: verify_operator('is null')
        True
        sage: verify_operator('not_an_operator')
        Traceback (most recent call last):
        ...
        TypeError: not_an_operator is not a legal operator.
    """
    binaries = ['=','<=','>=','like','<','>','<>','regexp']
    unaries = ['is null','is not null']
    if operator not in binaries and operator not in unaries:
        raise TypeError('%s is not a legal operator.'%operator)
    return True

def construct_skeleton(connection):
    """
    Constructs a database skeleton from the sql data.  The skeleton data
    structure is a triple indexed dictionary of the following format:

    skeleton -- a triple-indexed dictionary
        outer key - table name
            inner key - column name
                inner inner key - one of the following:
                    primary_key - boolean, whether column has been set as primary key
                    index - boolean, whether column has been set as index
                    sql - one of 'TEXT', 'BOOLEAN', 'INTEGER', 'REAL', or other
                            user defined type

    An example skeleton of a database with one table, that table with one column:
        {'table1':{'col1':{'primary_key':False, 'index':True, 'sql':'REAL'}}}

    EXAMPLE:
        sage: G = SQLDatabase(GraphDatabase().__dblocation__)
        sage: con = G.get_connection()
        sage: from sage.databases.database import construct_skeleton
        sage: construct_skeleton(con).keys()
        [u'aut_grp', u'degrees', u'spectrum', u'misc', u'graph_data']
    """
    skeleton = {}
    cur = connection.cursor()
    exe = cur.execute("select name from sqlite_master where type='table'")
    for table in exe.fetchall():
        skeleton[table[0]] = {}
        exe1 = cur.execute("pragma table_info(%s)"%table[0])
        for column in exe1.fetchall():
            skeleton[table[0]][column[1]] = {'sql':column[2], 'primary_key':(column[5]!=0), 'index':False}
        exe2 = cur.execute("pragma index_list(%s)"%table[0])
        for column in exe2.fetchall():
            if (column[1].find('sqlite') == -1):
                skeleton[table[0]][column[1]]['index'] = True
    return skeleton

def skel_to_col_attr_list(table_dict):
    """
    Returns a list of tuples representing the column attributes for a table.
    Each tuple represents one column and has values:
    (column_title <String>, sql_data_type <String>, primary_key <Boolean>)

    EXAMPLE:
        sage: from sage.databases.database import skel_to_col_attr_list, construct_skeleton
        sage: G = SQLDatabase(GraphDatabase().__dblocation__)
        sage: con = G.get_connection()
        sage: table_dict = construct_skeleton(con)['graph_data']
        sage: skel_to_col_attr_list(table_dict)
        [(u'perfect', u'BOOLEAN', False),
         (u'planar', u'BOOLEAN', False),
         (u'graph_id', u'INTEGER', False),
         (u'complement_graph6', u'TEXT', False),
         (u'num_edges', u'INTEGER', False),
         (u'num_cycles', u'INTEGER', False),
         (u'graph6', u'TEXT', False),
         (u'num_hamiltonian_cycles', u'INTEGER', False),
         (u'lovasz_number', u'REAL', False),
         (u'eulerian', u'BOOLEAN', False),
         (u'num_vertices', u'INTEGER', False)]

    """
    s = []
    for col in table_dict:
        s.append((col, table_dict[col]['sql'], table_dict[col]['primary_key']))
    return s

def new_table_set_col_attr(connection, table_name, table_skeleton):
    """
    Adds indices to the database wherever they are specified in the
    table_skeleton.  Note that an error will occur if adding an index
    on an already-indexed column.  We also remark that this is called
    automatically when creating a new table, so it is unlikely that a
    user will call this method directly.

    EXAMPLE:
        sage: S = SQLDatabase()
        sage: table_skeleton = {
        ...         'graph6':{'sql':'TEXT', 'index':True, 'primary_key':True},
        ...         'vertices':{'sql':'INTEGER'},
        ...         'edges':{'sql':'INTEGER'}
        ...         }
        sage: S.create_table('graph',table_skeleton)
        sage: diff_skeleton = {
        ...         'vertices':{'sql':'INTEGER', 'index':True},
        ...         'edges':{'sql':'INTEGER', 'index':True}
        ...         }
        sage: from sage.databases.database import new_table_set_col_attr, construct_skeleton
        sage: new_table_set_col_attr(S.get_connection(),'graph',diff_skeleton)
        sage: construct_skeleton(S.get_connection())
        {u'graph': {u'edges': {'index': True,
                               'primary_key': False,
                               'sql': u'INTEGER'},
                    u'graph6': {'index': True, 'primary_key': True, 'sql': u'TEXT'},
                    u'vertices': {'index': True,
                                  'primary_key': False,
                                  'sql': u'INTEGER'}}}
    """
    statement = ''
    for col in table_skeleton:
        if table_skeleton[col].has_key('index'):
            if table_skeleton[col]['index']:
                statement += 'CREATE INDEX %s ON %s (%s);\n'%(col, table_name, col)
            else:
                table_skeleton[col]['index'] = False
    if (statement != ''):
        connection.executescript(statement)


def _apply_format(func, x, html, id=None):
    """
    Applies a format to an output entry when showing query results.  For
    an example of use, see the GenericSQLQuery show function.

    INPUT:
        func -- the function that reformats the data
        x -- the data to reformat
        html -- whether the output display is in the notebook
        id -- if the function needs an additional identifier of the data

    EXAMPLES:
        sage: from sage.databases.database import _apply_format
        sage: def f(lower=''):
        ...     return lower.capitalize()
        sage: _apply_format(f,'capitalize me',False)
        'Capitalize me'
        sage: _apply_format(f,'capitalize me',True)
        '<td bgcolor=white align=center> Capitalize me </td>'
    """
    if id is not None:
        func_str = str(func(x,id))
    else:
        func_str = str(func(x))
    if html:
        return '<td bgcolor=white align=center> ' + func_str + ' </td>'
    else:
        return func_str

def _apply_plot(func, x, fig_num, figsize=[1,1], with_label=True,**kwds):
    """
    Saves a representative plot to be drawn in the database output.
    Allows passing of plot keywords.  Data x is assumed to be an
    identifier for the object's plotting function.  For an example of
    use, see the GenericSQLQuery show function.

    INPUT:
        func -- the function that plots the object
        x -- the data to pass to plotting function
        fig_num -- integer used to save plot data
        with_label -- whether or not to print data x with picture
        **kwds -- plotting keywords are allowed

    EXAMPLES:
        sage: from sage.databases.database import _apply_plot
        sage: def p(graph6):
        ...     g = Graph(str(graph6))
        ...     return g.plot(layout='circular', vertex_size=30, vertex_labels=False)
        sage: _apply_plot(p, 'C?', 0)
        '<td bgcolor=white align=center>C?<br><img src="cell://0.png"></td>'
    """
    p = func(x)
    p.save('%d.png'%fig_num, figsize=figsize)
    if with_label:
        return '<td bgcolor=white align=center>%s<br><img src="cell://%d.png"></td>'%(x,fig_num)
    else:
        return '<td bgcolor=white align=center><img src="cell://%d.png"></td>'%fig_num

class GenericSQLQuery(SageObject):

    def __init__(self, database, query_string, param_tuple=None):
        """
        A query for a SQLite database.

        INPUT:
            database -- a SQLDatabase or GenericSQLDatabase object
            query_string -- a string representing the SQL query
            param_tuple -- a tuple of strings - what to replace question marks in
                query_string with

        NOTE:
            This query class is generally intended for developers and more
            advanced users. It allows you to execute any query, and so may be
            considered unsafe...

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

        TUTORIAL:
        The SQLDatabase class is for interactively building databases intended for
        queries. This may sound redundant, but it is important. If you want a
        database intended for quick lookup of entries in very large tables, much
        like a hash table (such as a Python dictionary), a SQLDatabase may not be
        what you are looking for. The strength of SQLDatabases is in queries,
        searches through the database with complicated criteria.

        The class GenericSQLDatabase is for developers to provide a static
        database. The class does not support modification, and is meant to be a
        base class for specific classes of databases, such as the graph database.

        For example, we create a new database for storing isomorphism classes of
        simple graphs:
            sage: D = SQLDatabase()

        In order to generate representatives for the classes, we will import a
        function which generates all labeled graphs (noting that this is not the
        optimal way):
            sage: from sage.graphs.graph_isom import all_labeled_graphs

        We will need a table in the database in which to store the graphs, and we
        specify its structure with a Python dictionary, each of whose keys is the
        name of a column:
            sage: table_skeleton = {
            ... 'graph6':{'sql':'TEXT', 'index':True, 'primary_key':True},
            ... 'vertices':{'sql':'INTEGER'},
            ... 'edges':{'sql':'INTEGER'}
            ... }

        Then we create the table:
            sage: D.create_table('simon', table_skeleton)
            sage: D.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------

        Now that we have the table, we will begin to populate the table with
        rows. First, add the graph on zero vertices.
            sage: G = Graph()

            sage: D.add_row('simon',(0, G.graph6_string(), 0))

            sage: D.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------
            0                    ?                    0

        Next, add the graph on one vertex.
            sage: G.add_vertex()
            sage: D.add_row('simon',(0, G.graph6_string(), 1))
            sage: D.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------
            0                    ?                    0
            0                    @                    1

        Say we want a database of graphs on four or less vertices:
            sage: labels = {}
            sage: for i in range(2, 5):
            ...       labels[i] = []
            ...       for g in all_labeled_graphs(i):
            ...           g = g.canonical_label()
            ...           if g not in labels[i]:
            ...               labels[i].append(g)
            ...               D.add_row('simon', (g.size(), g.graph6_string(), g.order()))
            ...
            sage: D.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------
            0                    ?                    0
            0                    @                    1
            0                    A?                   2
            1                    A_                   2
            0                    B?                   3
            1                    BG                   3
            2                    BW                   3
            3                    Bw                   3
            0                    C?                   4
            1                    C@                   4
            2                    CB                   4
            3                    CF                   4
            3                    CJ                   4
            2                    CK                   4
            3                    CL                   4
            4                    CN                   4
            4                    C]                   4
            5                    C^                   4
            6                    C~                   4

        We can then query the database-- let's ask for all the graphs on four
        vertices with three edges. We do so by creating two queries and asking for
        rows that satisfy them both:
            sage: Q = SQLQuery(D, {'table_name':'simon', 'display_cols':['graph6'], 'expression':['vertices','=',4]})
            sage: Q2 = SQLQuery(D, {'table_name':'simon', 'display_cols':['graph6'], 'expression':['edges','=',3]})
            sage: Q = Q.intersect(Q2)
            sage: Q.run_query()
            [(u'CF', u'CF'), (u'CJ', u'CJ'), (u'CL', u'CL')]

        NOTE - The values of display_cols are always concatenated in intersections
        and unions.

        Of course, we can save the database to file:
            sage: replace_with_your_own_filepath = tmp_dir()
            sage: D.save(replace_with_your_own_filepath + 'simon.db')

        Now the database's hard link is to this file, and not the temporary db
        file. For example, let's say we open the same file with another class
        instance. We can load the file as an immutable database:
            sage: E = GenericSQLDatabase(replace_with_your_own_filepath + 'simon.db')
            sage: E.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------
            0                    ?                    0
            0                    @                    1
            0                    A?                   2
            1                    A_                   2
            0                    B?                   3
            1                    BG                   3
            2                    BW                   3
            3                    Bw                   3
            0                    C?                   4
            1                    C@                   4
            2                    CB                   4
            3                    CF                   4
            3                    CJ                   4
            2                    CK                   4
            3                    CL                   4
            4                    CN                   4
            4                    C]                   4
            5                    C^                   4
            6                    C~                   4
            sage: E.drop_table('simon')
            Traceback (most recent call last):
            ...
            AttributeError: 'GenericSQLDatabase' object has no attribute 'drop_table'
        """
        if not isinstance(database, GenericSQLDatabase):
            raise TypeError('%s is not a valid SQLDatabase'%database)
        self.__param_tuple__ = param_tuple
        self.__query_string__ = query_string
        self.__database__ = database

    def __repr__(self):
        """
        Overrides the print output to display useful info regarding the
        query.

        EXAMPLES:
            sage: G = GraphDatabase()
            sage: q = 'select graph_id,graph6,num_vertices,num_edges from graph_data where graph_id<=(?) and num_vertices=(?)'
            sage: param = (22,5)
            sage: Q = GenericSQLQuery(G,q,param)
            sage: print Q
            Query for sql database: ...graphs.db
            Query string: select graph_id,graph6,num_vertices,num_edges from graph_data where graph_id<=(?) and num_vertices=(?)
            Parameter tuple: (22, 5)

            sage: r = 'select graph6 from graph_data where num_vertices<=3'
            sage: R = GenericSQLQuery(G,r)
            sage: print R
            Query for sql database: ...graphs.db
            Query string: select graph6 from graph_data where num_vertices<=3
        """
        if not self.__query_string__: return 'Empty query on %s.'%self.__database__.__dblocation__
        s = "Query for sql database: "
        s += self.__database__.__dblocation__ + "\n"
        s += "Query string: "
        s += self.__query_string__ + "\n"
        if self.__param_tuple__:
            s += "Parameter tuple: "
            s += str(self.__param_tuple__) + "\n"
        return s

    def get_query_string(self):
        """
        Returns a copy of the query string.

        EXAMPLES:
            sage: G = GraphDatabase()
            sage: q = 'select graph_id,graph6,num_vertices,num_edges from graph_data where graph_id<=(?) and num_vertices=(?)'
            sage: param = (22,5)
            sage: Q = GenericSQLQuery(G,q,param)
            sage: Q.get_query_string()
            'select graph_id,graph6,num_vertices,num_edges from graph_data where graph_id<=(?) and num_vertices=(?)'
            sage: r = 'select graph6 from graph_data where num_vertices<=3'
            sage: R = GenericSQLQuery(G,r)
            sage: R.get_query_string()
            'select graph6 from graph_data where num_vertices<=3'
        """
        from copy import copy
        return copy(self.__query_string__)

    def copy(self):
        """
        Returns a copy of the query.

        EXAMPLES:
            sage: G = GraphDatabase()
            sage: q = 'select graph_id,graph6,num_vertices,num_edges from graph_data where graph_id<=(?) and num_vertices=(?)'
            sage: param = (22,5)
            sage: Q = GenericSQLQuery(G,q,param)
            sage: R = Q.copy()
            sage: R.get_query_string()
            'select graph_id,graph6,num_vertices,num_edges from graph_data where graph_id<=(?) and num_vertices=(?)'
            sage: R.show()
            graph_id             graph6               num_vertices         num_edges
            --------------------------------------------------------------------------------
            18                   D??                  5                    0
            19                   D?C                  5                    1
            20                   D?K                  5                    2
            21                   D@O                  5                    2
            22                   D?[                  5                    3
        """
        return GenericSQLQuery(self.__database__, self.__query_string__, self.__param_tuple__)

    def run_query(self):
        """
        Runs the query by executing the __query_string__.  Returns the results
        of the query in a list.

        EXAMPLES:
            sage: G = GraphDatabase()
            sage: q = 'select graph_id,graph6,num_vertices,num_edges from graph_data where graph_id<=(?) and num_vertices=(?)'
            sage: param = (22,5)
            sage: Q = GenericSQLQuery(G,q,param)
            sage: Q.run_query()
            [(18, u'D??', 5, 0),
             (19, u'D?C', 5, 1),
             (20, u'D?K', 5, 2),
             (21, u'D@O', 5, 2),
             (22, u'D?[', 5, 3)]

            sage: R = SQLQuery(G,{'table_name':'graph_data', 'display_cols':['graph6'], 'expression':['num_vertices','=',4]})
            sage: R.run_query()
            [(u'C?',),
             (u'C@',),
             (u'CB',),
             (u'CK',),
             (u'CF',),
             (u'CJ',),
             (u'CL',),
             (u'CN',),
             (u'C]',),
             (u'C^',),
             (u'C~',)]
        """
        try:
            cur = self.__database__.__connection__.cursor()
            if self.__param_tuple__ is not None:
                tup = []
                # make it a tuple of strings:
                for i in range (len(self.__param_tuple__)):
                    tup.append(str(self.__param_tuple__[i]))
                exe = cur.execute(self.__query_string__, tuple(tup))
            else:
                exe = cur.execute(self.__query_string__)
            lis = exe.fetchall()
            return lis
        except:
            raise RuntimeError('Failure to fetch query.')

    def show(self, max_field_size=20, format_cols={}, plot_cols={}, relabel_cols={}, id_col=None, **kwds):
        """
        Displays the result of the query in table format.

        INPUT:
            max_field_size -- how wide each field can be
            format_cols -- a dictionary that allows the user to specify the format
                            of a column's output by supplying a function.  The
                            format of the dictionary is:
                            {'column_name':(lambda x: format_function(x))}
            plot_cols -- a dictionary that allows the user to specify that a plot
                            should be drawn by the object generated by a data slice.
                            Note that plot kwds are permitted.  The dictionary format
                            is:
                            {'column_name':((lambda x: plot_function(x)), **kwds)}
            relabel_cols -- a dictionary to specify a relabeling of column headers.
                            The dictionary format is:
                            {'table_name':{'old_col_name':'new_col_name'}}
            id_col -- reference to a column that can be used as an object identifier
                            for each row

        EXAMPLE:
            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_data('simon',[(0,0),(1,1),(1,2)])
            sage: r = SQLQuery(DB, {'table_name':'simon', 'display_cols':['a1'], 'expression':['b2','<=', 6]})
            sage: p = SQLQuery(DB, {'table_name':'simon', 'display_cols':['b2'], 'expression':['b2','<=', 6]})
            sage: s = p.intersect(r)
            sage: s.show()
            b2                   a1
            ----------------------------------------
            0                    0
            1                    1
            2                    1

        sage: D = GraphDatabase()
        sage: from sage.graphs.graph_database import valid_kwds, data_to_degseq
        sage: relabel = {}
        sage: D = GraphDatabase()
        sage: from sage.graphs.graph_database import valid_kwds, data_to_degseq
        sage: relabel = {}
        sage: for col in valid_kwds:
        ...     relabel[col] = ' '.join([word.capitalize() for word in col.split('_')])

        sage: q = GraphQuery(display_cols=['graph6','degree_sequence'], num_vertices=4)
        sage: GenericSQLQuery.show(q, format_cols={'degree_sequence': (lambda x,y: data_to_degseq(x,y))}, relabel_cols=relabel,id_col='graph6')
        Graph6               Degree Sequence
        ----------------------------------------
        C?                   [0, 0, 0, 0]
        C@                   [0, 0, 1, 1]
        CB                   [0, 1, 1, 2]
        CK                   [1, 1, 1, 1]
        CF                   [1, 1, 1, 3]
        CJ                   [0, 2, 2, 2]
        CL                   [1, 1, 2, 2]
        CN                   [1, 2, 2, 3]
        C]                   [2, 2, 2, 2]
        C^                   [2, 2, 3, 3]
        C~                   [3, 3, 3, 3]
        """
        if not self.__query_string__:
            self.__database__.show()
            return

        try:
            cur = self.__database__.__connection__.cursor()
            if self.__param_tuple__ is not None:
                tup = []
                # make it a tuple of strings:
                for i in range (len(self.__param_tuple__)):
                    tup.append(str(self.__param_tuple__[i]))
                cur.execute(self.__query_string__, tuple(tup))
            else:
                cur.execute(self.__query_string__)
        except:
            raise RuntimeError('Failure to fetch query.')

        col_titles = [des[0] for des in cur.description]
        fcol_index = []
        fcol_titles = []
        pcol_index = []
        pcol_titles = []

        if id_col:
            id_col_index = col_titles.index(id_col)
        else:
            id_col_index = None

        if format_cols:
            for col in format_cols:
                fcol_index.append(col_titles.index(col))
                fcol_titles.append(col)

        if plot_cols:
            for col in plot_cols:
                pcol_index.append(col_titles.index(col))
                pcol_titles.append(col)

        if relabel_cols:
            for i in range(len(col_titles)):
                try:
                    col_titles[i] = relabel_cols[col_titles[i]]
                except KeyError:
                    continue

        from sage.server.support import EMBEDDED_MODE
        if EMBEDDED_MODE:
            # Notebook Version
            print '<html><!--notruncate--><table bgcolor=lightgrey cellpadding=0><tr>'
            for col in col_titles:
                print '<td bgcolor=white align=center> ' + col + ' </td>'
            print '</tr>'
            p = 0
            for row in cur:
                f = 0
                print '<tr>'
                for index in range(len(col_titles)):
                    if index in fcol_index:
                        if id_col_index is not None:
                            print _apply_format(format_cols[fcol_titles[f]], row[index], True, row[id_col_index])
                        else:
                            print _apply_format(format_cols[fcol_titles[f]], row[index], True)
                        f += 1
                    elif index in pcol_index:
                        print _apply_plot(plot_cols[pcol_titles[p%len(pcol_titles)]], row[index], p)
                        p += 1
                    else:
                        print '<td bgcolor=white align=center> ' + str(row[index]) + ' </td>'
                print '</tr>'
            print '</table></html>'

        else:
            # Command Prompt Version
            for col in col_titles:
                print col.ljust(max_field_size),
            print # new line
            print '-' * max_field_size * len(col_titles)


            for row in cur:
                f = 0
                for index in range(len(col_titles)):
                    if index in fcol_index:
                        if id_col_index is not None:
                            field_val = _apply_format(format_cols[fcol_titles[f]], \
                                                    row[index], False, row[id_col_index])
                        else:
                            field_val = _apply_format(format_cols[fcol_titles[f]], \
                                                    row[index], False)
                        f += 1
                    elif index in pcol_index:
                        raise NotImplementedError, 'Cannot display plot on command line.'
                    else:
                        field_val = str(row[index])
                    print field_val.ljust(max_field_size).ljust(max_field_size),
                print # new line

class SQLQuery(GenericSQLQuery):
    def __init__(self, database, query_dict={}):
        """
        A query for a SQLite database.

        INPUT:
            database -- a SQLDatabase or GenericSQLDatabase object
            query_dict -- a dictionary specifying the query itself. The format is:

        {'table_name': 'tblname', 'display_cols': ['col1', 'col2', 'col3'], 'expression':[col, operator, value]}

        NOTE:
            Every SQL type we are using is ultimately represented as a string, so
            if you wish to save actual strings to a database, you actually need to
            do something like: '"value"'.

        TUTORIAL:
        The SQLDatabase class is for interactively building databases intended for
        queries. This may sound redundant, but it is important. If you want a
        database intended for quick lookup of entries in very large tables, much
        like a hash table (such as a Python dictionary), a SQLDatabase may not be
        what you are looking for. The strength of SQLDatabases is in queries,
        searches through the database with complicated criteria.

        The class GenericSQLDatabase is for developers to provide a static
        database. The class does not support modification, and is meant to be a
        base class for specific classes of databases, such as the graph database.

        For example, we create a new database for storing isomorphism classes of
        simple graphs:
            sage: D = SQLDatabase()

        In order to generate representatives for the classes, we will import a
        function which generates all labeled graphs (noting that this is not the
        optimal way):
            sage: from sage.graphs.graph_isom import all_labeled_graphs

        We will need a table in the database in which to store the graphs, and we
        specify its structure with a Python dictionary, each of whose keys is the
        name of a column:
            sage: table_skeleton = {
            ... 'graph6':{'sql':'TEXT', 'index':True, 'primary_key':True},
            ... 'vertices':{'sql':'INTEGER'},
            ... 'edges':{'sql':'INTEGER'}
            ... }

        Then we create the table:
            sage: D.create_table('simon', table_skeleton)
            sage: D.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------

        Now that we have the table, we will begin to populate the table with
        rows. First, add the graph on zero vertices.
            sage: G = Graph()

            sage: D.add_row('simon',(0, G.graph6_string(), 0))

            sage: D.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------
            0                    ?                    0

        Next, add the graph on one vertex.
            sage: G.add_vertex()
            sage: D.add_row('simon',(0, G.graph6_string(), 1))
            sage: D.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------
            0                    ?                    0
            0                    @                    1

        Say we want a database of graphs on four or less vertices:
            sage: labels = {}
            sage: for i in range(2, 5):
            ...       labels[i] = []
            ...       for g in all_labeled_graphs(i):
            ...           g = g.canonical_label()
            ...           if g not in labels[i]:
            ...               labels[i].append(g)
            ...               D.add_row('simon', (g.size(), g.graph6_string(), g.order()))
            ...
            sage: D.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------
            0                    ?                    0
            0                    @                    1
            0                    A?                   2
            1                    A_                   2
            0                    B?                   3
            1                    BG                   3
            2                    BW                   3
            3                    Bw                   3
            0                    C?                   4
            1                    C@                   4
            2                    CB                   4
            3                    CF                   4
            3                    CJ                   4
            2                    CK                   4
            3                    CL                   4
            4                    CN                   4
            4                    C]                   4
            5                    C^                   4
            6                    C~                   4

        We can then query the database-- let's ask for all the graphs on four
        vertices with three edges. We do so by creating two queries and asking for
        rows that satisfy them both:
            sage: Q = SQLQuery(D, {'table_name':'simon', 'display_cols':['graph6'], 'expression':['vertices','=',4]})
            sage: Q2 = SQLQuery(D, {'table_name':'simon', 'display_cols':['graph6'], 'expression':['edges','=',3]})
            sage: Q = Q.intersect(Q2)
            sage: Q.run_query()
            [(u'CF', u'CF'), (u'CJ', u'CJ'), (u'CL', u'CL')]

        NOTE - The values of display_cols are always concatenated in intersections
        and unions.

        Of course, we can save the database to file:
            sage: replace_with_your_own_filepath = tmp_dir()
            sage: D.save(replace_with_your_own_filepath + 'simon.db')

        Now the database's hard link is to this file, and not the temporary db
        file. For example, let's say we open the same file with another class
        instance. We can load the file as an immutable database:
            sage: E = GenericSQLDatabase(replace_with_your_own_filepath + 'simon.db')
            sage: E.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------
            0                    ?                    0
            0                    @                    1
            0                    A?                   2
            1                    A_                   2
            0                    B?                   3
            1                    BG                   3
            2                    BW                   3
            3                    Bw                   3
            0                    C?                   4
            1                    C@                   4
            2                    CB                   4
            3                    CF                   4
            3                    CJ                   4
            2                    CK                   4
            3                    CL                   4
            4                    CN                   4
            4                    C]                   4
            5                    C^                   4
            6                    C~                   4
            sage: E.drop_table('simon')
            Traceback (most recent call last):
            ...
            AttributeError: 'GenericSQLDatabase' object has no attribute 'drop_table'

        """
        self.__database__ = database
        if not query_dict:
            self.__query_dict__ = {}
        else:
            if not database.__skeleton__.has_key(query_dict['table_name']):
                raise ValueError("Database has no table %s."%query_dict['table_name'])
            if query_dict['display_cols'] is not None:
                for column in query_dict['display_cols']:
                    if not database.__skeleton__[query_dict['table_name']].has_key(column):
                        raise ValueError("Table has no column %s."%column)
            if not database.__skeleton__[query_dict['table_name']].has_key(query_dict['expression'][0]):
                raise ValueError("Table has no column %s."%query_dict['expression'][0])
            self.__query_dict__ = query_dict

        if not isinstance(database, GenericSQLDatabase):
            raise TypeError('%s is not a valid SQLDatabase'%database)

        self.__database__ = database

        # confirm operator:
        if self.__query_dict__:
            verify_operator(query_dict['expression'][1])

        # make tuple:
        if self.__query_dict__:
            self.__param_tuple__ = (self.__query_dict__['expression'][2],)
        else:
            self.__param_tuple__ = tuple()

        # make query string:
        if self.__query_dict__:
            # display cols:
            if self.__query_dict__['display_cols'] is not None:
                for i in range(len(self.__query_dict__['display_cols'])):
                    self.__query_dict__['display_cols'][i] = self.__query_dict__['table_name'] + '.'+ self.__query_dict__['display_cols'][i]
                self.__query_string__ = 'SELECT ' + ', '.join(self.__query_dict__['display_cols']) + \
                                    ' FROM ' + self.__query_dict__['table_name'] + \
                                    ' WHERE ' + self.__query_dict__['table_name'] + '.' + \
                                    self.__query_dict__['expression'][0] + ' ' + \
                                    self.__query_dict__['expression'][1] + ' ?'
            else:
                self.__query_string__ = 'SELECT ' + ', ' + \
                                    ' FROM ' + self.__query_dict__['table_name'] + \
                                    ' WHERE ' + self.__query_dict__['table_name'] + '.' + \
                                    self.__query_dict__['expression'][0] + ' ' + \
                                    self.__query_dict__['expression'][1] + ' ?'
        else:
            self.__query_string__ = ''

    def copy(self):
        """
        Returns a copy of itself.

        EXAMPLES:
            sage: G = GraphDatabase()
            sage: Q = SQLQuery(G, {'table_name':'graph_data', 'display_cols':['graph_id','graph6','num_vertices'],'expression':['num_edges','<',3]})
            sage: Q.show()
            graph_id             graph6               num_vertices
            ------------------------------------------------------------
            0                    @                    1
            1                    A?                   2
            3                    B?                   3
            7                    C?                   4
            18                   D??                  5
            52                   E???                 6
            208                  F????                7
            2                    A_                   2
            4                    BG                   3
            8                    C@                   4
            19                   D?C                  5
            53                   E??G                 6
            209                  F???G                7
            5                    BW                   3
            9                    CB                   4
            10                   CK                   4
            20                   D?K                  5
            21                   D@O                  5
            54                   E??W                 6
            55                   E?C_                 6
            210                  F???W                7
            211                  F??G_                7
            sage: R = Q.copy()
            sage: print R
            Query for sql database: ...graphs.db
            Query string: SELECT graph_data.graph_id, graph_data.graph6, graph_data.num_vertices FROM graph_data WHERE graph_data.num_edges < ?
            Parameter tuple: (3,)

            sage: R.show()
            graph_id             graph6               num_vertices
            ------------------------------------------------------------
            0                    @                    1
            1                    A?                   2
            3                    B?                   3
            7                    C?                   4
            18                   D??                  5
            52                   E???                 6
            208                  F????                7
            2                    A_                   2
            4                    BG                   3
            8                    C@                   4
            19                   D?C                  5
            53                   E??G                 6
            209                  F???G                7
            5                    BW                   3
            9                    CB                   4
            10                   CK                   4
            20                   D?K                  5
            21                   D@O                  5
            54                   E??W                 6
            55                   E?C_                 6
            210                  F???W                7
            211                  F??G_                7
        """
        d = SQLQuery(self.__database__)
        d.__query_dict__ = self.__query_dict__
        d.__query_string__ = self.__query_string__
        d.__param_tuple__ = self.__param_tuple__
        return d

    def intersect(self, other, join_table=None, join_dict=None, in_place=False):
        """
        Returns a new SQLQuery that is the intersection of self and other.
        join_table and join_dict can be None iff the two queries only search
        one table in the database.  All display columns will be concatenated in
        order: self display cols + other display cols.

        INPUT:
            other -- the SQLQuery to intersect with
            join_table -- base table to join on (This table should have at least
                one column in each table to join on).
            join_dict -- a dictionary that represents the join structure for the
                new query.  (Must include a mapping for all tables, including
                those previously joined in either query).  Structure is given
                by:
                    {'join_table1': ('corr_base_col1', 'col1'), 'join_table2': ('corr_base_col2', 'col2')}
                where join_table1 is to be joined with join_table on
                    join_table.corr_base_col1 = join_table1.col1

        EXAMPLES:
            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.create_table('lucy',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_data('simon', [(0,5),(1,4)])
            sage: DB.add_data('lucy', [(1,1),(1,4)])
            sage: q = SQLQuery(DB, {'table_name':'lucy', 'display_cols':['b2'], 'expression':['a1','=',1]})
            sage: r = SQLQuery(DB, {'table_name':'simon', 'display_cols':['a1'], 'expression':['b2','<=', 6]})
            sage: s = q.intersect(r, 'simon', {'lucy':('a1','a1')})
            sage: s.run_query()
            [(1, 1), (4, 1)]
            sage: s = q.intersect(r)
            Traceback (most recent call last):
            ...
            ValueError: Input queries query different tables but join parameters are NoneType

        """
        if self.__database__.__dblocation__ != other.__database__.__dblocation__:
            raise TypeError('Queries %s and %s must be attached to the same database.'%(self, other))

        if in_place:

            if not self.__query_string__:
                self.__query_string__ = other.__query_string__
                self.__param_tuple__ = other.__param_tuple__
                return
            if not other.__query_string__: return
            if join_table is None or join_dict is None:
                pattern = ' JOIN '
                if re.search(pattern,self.__query_string__) or re.search(pattern,other.__query_string__):
                    raise TypeError('Input queries have joins but join parameters are NoneType.')
                s = ((self.__query_string__).upper()).split('FROM ')
                o = ((other.__query_string__).upper()).split('FROM ')
                s = s[1].split(' WHERE ')
                o = o[1].split(' WHERE ')
                if s[0] != o[0]:
                    raise ValueError('Input queries query different tables but join parameters are NoneType')

            # inner join clause
            if join_dict is not None:
                joins = join_table
                for table in join_dict:
                    joins += ' INNER JOIN ' + table + ' ON ' + join_table + \
                             '.' + join_dict[table][0] + '=' + table + '.' + \
                             join_dict[table][1] + ' '
                self.__query_string__ = re.sub(' FROM .* WHERE ', ' FROM ' + joins + 'WHERE ', self.__query_string__)

            # concatenate display cols
            disp = self.__query_string__.split(' FROM')
            disp[0] += ',' + other.__query_string__.split(' FROM')[0].split('SELECT ')[1]+' FROM'
            new_query = ''.join(disp)

            # concatenate where clause
            new_query = re.sub(' WHERE ',' WHERE ( ',new_query)
            new_query += re.sub('^.* WHERE ',' ) AND ( ',other.__query_string__)
            self.__query_string__ = new_query + ' )'

            if self.__param_tuple__ is None or not self.__param_tuple__:
                self.__param_tuple__ = other.__param_tuple__
            elif other.__param_tuple__ is None or not other.__param_tuple__:
                self.__param_tuple__ = self.__param_tuple__
            else:
                self.__param_tuple__ = self.__param_tuple__ + other.__param_tuple__

        else:
            if not self.__query_string__: return other.copy()
            if not other.__query_string__: return self.copy()
            if join_table is None or join_dict is None:
                pattern = ' JOIN '
                if re.search(pattern,self.__query_string__) or re.search(pattern,other.__query_string__):
                    raise TypeError('Input queries have joins but join parameters are NoneType.')
                s = ((self.__query_string__).upper()).split('FROM ')
                o = ((other.__query_string__).upper()).split('FROM ')
                s = s[1].split(' WHERE ')
                o = o[1].split(' WHERE ')
                if s[0] != o[0]:
                    raise ValueError('Input queries query different tables but join parameters are NoneType')

            q = self.copy()

            # inner join clause
            if join_dict is not None:
                joins = join_table
                for table in join_dict:
                    joins += ' INNER JOIN ' + table + ' ON ' + join_table + \
                             '.' + join_dict[table][0] + '=' + table + '.' + \
                             join_dict[table][1] + ' '
                q.__query_string__ = re.sub(' FROM .* WHERE ', ' FROM ' + joins + 'WHERE ', self.__query_string__)

            # concatenate display cols
            disp = q.__query_string__.split(' FROM')
            disp[0] += ',' + other.__query_string__.split(' FROM')[0].split('SELECT ')[1]+' FROM'
            new_query = ''.join(disp)

            # concatenate where clause
            new_query = re.sub(' WHERE ',' WHERE ( ',new_query)
            new_query += re.sub('^.* WHERE ',' ) AND ( ',other.__query_string__)
            q.__query_string__ = new_query + ' )'

            if self.__param_tuple__ is None or not self.__param_tuple__:
                q.__param_tuple__ = other.__param_tuple__
            elif other.__param_tuple__ is None or not other.__param_tuple__:
                q.__param_tuple__ = self.__param_tuple__
            else:
                q.__param_tuple__ = self.__param_tuple__ + other.__param_tuple__

            return q

    def union(self, other, join_table=None, join_dict=None, in_place=False):
        """
        Returns a new SQLQuery that is the union of self and other.
        join_table and join_dict can be None iff the two queries only search
        one table in the database.  All display columns will be concatenated in
        order: self display cols + other display cols.

        INPUT:
            other -- the SQLQuery to union with
            join_table -- base table to join on (This table should have at least
                one column in each table to join on).
            join_dict -- a dictionary that represents the join structure for the
                new query.  (Must include a mapping for all tables, including
                those previously joined in either query).  Structure is given
                by:
                    {'join_table1': ('corr_base_col1', 'col1'), 'join_table2': ('corr_base_col2', 'col2')}
                where join_table1 is to be joined with join_table on
                    join_table.corr_base_col1 = join_table1.col1

        EXAMPLES:
            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.create_table('lucy',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_data('simon', [(0,5),(1,4)])
            sage: DB.add_data('lucy', [(1,1),(1,4)])
            sage: q = SQLQuery(DB, {'table_name':'lucy', 'display_cols':['b2'], 'expression':['a1','=',1]})
            sage: r = SQLQuery(DB, {'table_name':'simon', 'display_cols':['a1'], 'expression':['b2','<=', 6]})
            sage: s = q.union(r, 'simon', {'lucy':('a1','a1')})
            sage: s.__query_string__
            'SELECT lucy.b2,simon.a1 FROM simon INNER JOIN lucy ON simon.a1=lucy.a1 WHERE ( lucy.a1 = ? ) OR ( simon.b2 <= ? )'

            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.create_table('lucy',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_data('simon', [(0,5),(1,4)])
            sage: DB.add_data('lucy', [(1,1),(1,4)])
            sage: q = SQLQuery(DB, {'table_name':'lucy', 'display_cols':['b2'], 'expression':['a1','=',1]})
            sage: r = SQLQuery(DB, {'table_name':'simon', 'display_cols':['a1'], 'expression':['b2','<=', 6]})
            sage: s = q.union(r, 'simon', {'lucy':('a1','a1')})
            sage: s.run_query()
            [(1, 1), (4, 1)]
            sage: s.show()
            b2                   a1
            ----------------------------------------
            1                    1
            4                    1
            sage: s = q.union(r)
            Traceback (most recent call last):
            ...
            ValueError: Input queries query different tables but join parameters are NoneType


        """
        if self.__database__.__dblocation__ != other.__database__.__dblocation__:
            raise TypeError('Queries %s and %s must be attached to the same database.'%(self, other))

        if in_place:
            if not self.__query_string__: return
            if not other.__query_string__: return
            if join_table is None or join_dict is None:
                pattern = ' JOIN '
                if re.search(pattern,self.__query_string__) or re.search(pattern,other.__query_string__):
                    raise TypeError('Input queries have joins but join parameters are NoneType.')
                s = ((self.__query_string__).upper()).split('FROM ')
                o = ((other.__query_string__).upper()).split('FROM ')
                s = s[1].split(' WHERE ')
                o = o[1].split(' WHERE ')
                if s[0] != o[0]:
                    raise ValueError('Input queries query different tables but join parameters are NoneType')

            # inner join clause
            if join_dict is not None:
                joins = join_table
                for table in join_dict:
                    joins += ' INNER JOIN ' + table + ' ON ' + join_table + \
                             '.' + join_dict[table][0] + '=' + table + '.' + \
                             join_dict[table][1] + ' '
                self.__query_string__ = re.sub(' FROM .* WHERE ', ' FROM ' + joins + 'WHERE ', self.__query_string__)

            # concatenate display cols
            disp = self.__query_string__.split(' FROM')
            disp[0] += ',' + other.__query_string__.split(' FROM')[0].split('SELECT ')[1]+' FROM'
            new_query = ''.join(disp)

            # concatenate where clause
            new_query = re.sub(' WHERE ',' WHERE ( ',new_query)
            new_query += re.sub('^.* WHERE ',' ) OR ( ',other.__query_string__)
            self.__query_string__ = new_query + ' )'

            self.__param_tuple__ = self.__param_tuple__ + other.__param_tuple__
        else:
            q = self.copy()
            if not self.__query_string__: return q
            if not other.__query_string__: return other.copy()
            if join_table is None or join_dict is None:
                pattern = ' JOIN '
                if re.search(pattern,self.__query_string__) or re.search(pattern,other.__query_string__):
                    raise TypeError('Input queries have joins but join parameters are NoneType.')
                s = ((self.__query_string__).upper()).split('FROM ')
                o = ((other.__query_string__).upper()).split('FROM ')
                s = s[1].split(' WHERE ')
                o = o[1].split(' WHERE ')
                if s[0] != o[0]:
                    raise ValueError('Input queries query different tables but join parameters are NoneType')

            # inner join clause
            if join_dict is not None:
                joins = join_table
                for table in join_dict:
                    joins += ' INNER JOIN ' + table + ' ON ' + join_table + \
                             '.' + join_dict[table][0] + '=' + table + '.' + \
                             join_dict[table][1] + ' '
                q.__query_string__ = re.sub(' FROM .* WHERE ', ' FROM ' + joins + 'WHERE ', self.__query_string__)

            # concatenate display cols
            disp = q.__query_string__.split(' FROM')
            disp[0] += ',' + other.__query_string__.split(' FROM')[0].split('SELECT ')[1]+' FROM'
            new_query = ''.join(disp)

            # concatenate where clause
            new_query = re.sub(' WHERE ',' WHERE ( ',new_query)
            new_query += re.sub('^.* WHERE ',' ) OR ( ',other.__query_string__)
            q.__query_string__ = new_query + ' )'

            q.__param_tuple__ = self.__param_tuple__ + other.__param_tuple__

            return q

    def complement(self, in_place=False):
        """
        Returns a new SQLQuery that is the complement of self.

        EXAMPLE:
            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_data('simon', [(0,5),(1,4)])
            sage: DB.show('simon')
            a1                   b2
            ----------------------------------------
            0                    5
            1                    4
            sage: r = SQLQuery(DB, {'table_name':'simon', 'display_cols':['a1'], 'expression':['b2','<=', 6]})
            sage: r.show()
            a1
            --------------------
            0
            1
            sage: s = r.complement()
            sage: s.show()
            a1
            --------------------

        """
        if in_place:
            if self.__query_string__:
                self.__query_string__ = re.sub(' WHERE ',' WHERE NOT ( ',self.__query_string__) + ' )'
            else: return
        else:
            if not self.__query_string__: return self.copy()
            q = SQLQuery(self.__database__)
            q.__query_string__ = re.sub(' WHERE ',' WHERE NOT ( ',self.__query_string__)
            q.__query_string__ += ' )'
            q.__param_tuple__ = self.__param_tuple__
            return q

class GenericSQLDatabase(SageObject):

    def __init__(self, filename):
        """
        *Immutable* Database class.

        INPUT:
            filename -- where to load the database from

        TUTORIAL:
        The SQLDatabase class is for interactively building databases intended for
        queries. This may sound redundant, but it is important. If you want a
        database intended for quick lookup of entries in very large tables, much
        like a hash table (such as a Python dictionary), a SQLDatabase may not be
        what you are looking for. The strength of SQLDatabases is in queries,
        searches through the database with complicated criteria.

        The class GenericSQLDatabase is for developers to provide a static
        database. The class does not support modification, and is meant to be a
        base class for specific classes of databases, such as the graph database.

        For example, we create a new database for storing isomorphism classes of
        simple graphs:
            sage: D = SQLDatabase()

        In order to generate representatives for the classes, we will import a
        function which generates all labeled graphs (noting that this is not the
        optimal way):
            sage: from sage.graphs.graph_isom import all_labeled_graphs

        We will need a table in the database in which to store the graphs, and we
        specify its structure with a Python dictionary, each of whose keys is the
        name of a column:
            sage: table_skeleton = {
            ... 'graph6':{'sql':'TEXT', 'index':True, 'primary_key':True},
            ... 'vertices':{'sql':'INTEGER'},
            ... 'edges':{'sql':'INTEGER'}
            ... }

        Then we create the table:
            sage: D.create_table('simon', table_skeleton)
            sage: D.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------

        Now that we have the table, we will begin to populate the table with
        rows. First, add the graph on zero vertices.
            sage: G = Graph()

            sage: D.add_row('simon',(0, G.graph6_string(), 0))

            sage: D.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------
            0                    ?                    0

        Next, add the graph on one vertex.
            sage: G.add_vertex()
            sage: D.add_row('simon',(0, G.graph6_string(), 1))
            sage: D.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------
            0                    ?                    0
            0                    @                    1

        Say we want a database of graphs on four or less vertices:
            sage: labels = {}
            sage: for i in range(2, 5):
            ...       labels[i] = []
            ...       for g in all_labeled_graphs(i):
            ...           g = g.canonical_label()
            ...           if g not in labels[i]:
            ...               labels[i].append(g)
            ...               D.add_row('simon', (g.size(), g.graph6_string(), g.order()))
            ...
            sage: D.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------
            0                    ?                    0
            0                    @                    1
            0                    A?                   2
            1                    A_                   2
            0                    B?                   3
            1                    BG                   3
            2                    BW                   3
            3                    Bw                   3
            0                    C?                   4
            1                    C@                   4
            2                    CB                   4
            3                    CF                   4
            3                    CJ                   4
            2                    CK                   4
            3                    CL                   4
            4                    CN                   4
            4                    C]                   4
            5                    C^                   4
            6                    C~                   4

        We can then query the database-- let's ask for all the graphs on four
        vertices with three edges. We do so by creating two queries and asking for
        rows that satisfy them both:
            sage: Q = SQLQuery(D, {'table_name':'simon', 'display_cols':['graph6'], 'expression':['vertices','=',4]})
            sage: Q2 = SQLQuery(D, {'table_name':'simon', 'display_cols':['graph6'], 'expression':['edges','=',3]})
            sage: Q = Q.intersect(Q2)
            sage: Q.run_query()
            [(u'CF', u'CF'), (u'CJ', u'CJ'), (u'CL', u'CL')]

        NOTE - The values of display_cols are always concatenated in intersections
        and unions.

        Of course, we can save the database to file:
            sage: replace_with_your_own_filepath = tmp_dir()
            sage: D.save(replace_with_your_own_filepath + 'simon.db')

        Now the database's hard link is to this file, and not the temporary db
        file. For example, let's say we open the same file with another class
        instance. We can load the file as an immutable database:
            sage: E = GenericSQLDatabase(replace_with_your_own_filepath + 'simon.db')
            sage: E.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------
            0                    ?                    0
            0                    @                    1
            0                    A?                   2
            1                    A_                   2
            0                    B?                   3
            1                    BG                   3
            2                    BW                   3
            3                    Bw                   3
            0                    C?                   4
            1                    C@                   4
            2                    CB                   4
            3                    CF                   4
            3                    CJ                   4
            2                    CK                   4
            3                    CL                   4
            4                    CN                   4
            4                    C]                   4
            5                    C^                   4
            6                    C~                   4
            sage: E.drop_table('simon')
            Traceback (most recent call last):
            ...
            AttributeError: 'GenericSQLDatabase' object has no attribute 'drop_table'


        """
        if (filename[-3:] != '.db'):
            raise ValueError('Please enter a valid database path (file name %s does not end in .db).'%filename)
        self.__dblocation__ = filename
        self.__connection__ = sqlite.connect(self.__dblocation__, check_same_thread=False)
        # check_same_thread = False:
        # this is to avoid the multiple thread problem with dsage:
        # pysqlite does not trust multiple threads for the same connection
        self.__connection__.create_function("regexp", 2, regexp)

        self.__skeleton__ = construct_skeleton(self.__connection__)

    def __repr__(self):
        """
        Overrides the print output to display useful info regarding the
        database.

        EXAMPLE:
            sage: replace_with_filepath = tmp_dir() + 'test.db'
            sage: SD = SQLDatabase(replace_with_filepath)
            sage: SD.create_table('simon', {'n':{'sql':'INTEGER', 'index':True}})
            sage: print SD
            table simon:
               column n: index: True; primary_key: False; sql: INTEGER;
        """
        s = ''
        for table in self.__skeleton__:
            s += 'table ' + table + ':\n'
            for column in self.__skeleton__[table]:
                s += '   column ' + column + ': '
                for data in self.__skeleton__[table][column]:
                    s += data + ': ' + str(self.__skeleton__[table][column][data]) + '; '
                s += '\n'
        return s

    def copy(self):
        """
        Returns an instance of SQLDatabase that points to a copy database,
        and allows modification.

        EXAMPLES:
            sage: DB = SQLDatabase()
            sage: DB.create_table('lucy',{'id':{'sql':'INTEGER', 'primary_key':True, 'index':True},'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_rows('lucy', [(0,1,1),(1,1,4),(2,0,7),(3,1,384),(4,1,978932)],['id','a1','b2'])
            sage: d = DB.copy()

            sage: d.show('lucy')
            a1                   id                   b2
            ------------------------------------------------------------
            1                    0                    1
            1                    1                    4
            0                    2                    7
            1                    3                    384
            1                    4                    978932

            sage: DB.show('lucy')
            a1                   id                   b2
            ------------------------------------------------------------
            1                    0                    1
            1                    1                    4
            0                    2                    7
            1                    3                    384
            1                    4                    978932

            sage: Q = SQLQuery(DB, {'table_name':'lucy', 'display_cols':['id','a1','b2'], 'expression':['id','>=',3]})
            sage: DB.delete_rows(Q)
            sage: DB.show('lucy')
            a1                   id                   b2
            ------------------------------------------------------------
            1                    0                    1
            1                    1                    4
            0                    2                    7

            sage: d.show('lucy')
            a1                   id                   b2
            ------------------------------------------------------------
            1                    0                    1
            1                    1                    4
            0                    2                    7
            1                    3                    384
            1                    4                    978932

        """
        from copy import copy
        # copy .db file
        new_loc = tmp_filename() + '.db'
        os.system('cp '+ self.__dblocation__ + ' ' + new_loc)
        D = SQLDatabase(filename=new_loc)
        for table in D.__skeleton__:
            # Get an ordered list:
            cur_list = skel_to_col_attr_list(D.__skeleton__[table])

            new = ''
            for col in cur_list:
                new += str(col[0]) +', '
            new = new.rstrip(', ')

            data = ((self.__connection__).execute('SELECT %s from %s'%(new,table))).fetchall()
            new = new.split(', ')

            # Fill data in new table
            D.add_rows(table_name=table,rows=data,entry_order=new)
        return D

    def save(self, filename):
        """
        Save the database to the specified location.

        EXAMPLE:
            sage: MonicPolys = SQLDatabase()
            sage: MonicPolys.create_table('simon', {'n':{'sql':'INTEGER', 'index':True}})
            sage: for n in range(20):
            ...     MonicPolys.add_row('simon', (n,))
            sage: tmp = tmp_dir() # replace with your own file path
            sage: MonicPolys.save(tmp+'sage.db')
            sage: N = GenericSQLDatabase(tmp+'sage.db')
            sage: N.show('simon')
            n
            --------------------
            0
            1
            2
            3
            4
            5
            6
            7
            8
            9
            10
            11
            12
            13
            14
            15
            16
            17
            18
            19
        """
        try:
            self.__connection__.execute('commit')
        except:
            # Not sure why this throws an exception - but without it,
            #       the changes are not committed so it is necessary.
            pass
        os.system('cp ' + self.__dblocation__ + ' ' + filename)

    def get_skeleton(self, check=False):
        """
        Returns a dictionary representing the hierarchical structure of the
        database, in the following format.

            skeleton -- a triple-indexed dictionary
                outer key - table name
                    inner key - column name
                        inner inner key - one of the following:
                primary_key - boolean, whether column has been set as primary key
                index - boolean, whether column has been set as index
                sql - one of 'TEXT', 'BOOLEAN', 'INTEGER', 'REAL', or other
                    user defined type

        For example,
        {'table1':{'col1':{'primary_key':False, 'index':True, 'sql':'REAL'}}}

        INPUT:
            check -- if True, checks to make sure the database's actual structure
            matches the skeleton on record.

        EXAMPLES:
            sage: GDB = GraphDatabase()
            sage: GDB.get_skeleton()             # slightly random output
            {u'aut_grp': {u'aut_grp_size': {'index': True,
                                'primary_key': False,
                                'sql': u'INTEGER'},
                                ...
               u'spectrum': {'index': False,
                             'primary_key': False,
                             'sql': u'TEXT'}}}

        """
        if not self.__skeleton__:
            self.__skeleton__ = construct_skeleton(self.__connection__)
        elif check:
            d = construct_skeleton(self.__connection__)
            if d == self.__skeleton__:
                return d
            else:
                raise RuntimeError("Skeleton structure is out of whack!")
        return self.__skeleton__

    def query(self, query_dict={}):
        """
        Creates a SQLQuery on this database.  For full class details,
        type SQLQuery? and press shift+enter.

        EXAMPLE:
            sage: D = SQLDatabase()
            sage: D.create_table('simon', {'wolf':{'sql':'BOOLEAN'}, 'tag':{'sql':'INTEGER'}})
            sage: q = D.query({'table_name':'simon', 'display_cols':['tag'], 'expression':['wolf','=',1]})
            sage: q.get_query_string()
            'SELECT simon.tag FROM simon WHERE simon.wolf = ?'
            sage: q.__param_tuple__
            (1,)
        """
        return SQLQuery(self, query_dict)

    def generic_query(self, query_string, param_tuple=None):
        """
        Creates a GenericSQLQuery on this database.  For full class details,
        type GenericSQLQuery? and press shift+enter.

        EXAMPLE:
            sage: D = SQLDatabase()
            sage: D.create_table('simon', {'wolf':{'sql':'BOOLEAN'}, 'tag':{'sql':'INTEGER'}})
            sage: q = D.generic_query(query_string='select tag from simon where wolf=?',param_tuple=(1,))
            sage: q.get_query_string()
            'select tag from simon where wolf=?'
            sage: q.__param_tuple__
            (1,)
        """
        return GenericSQLQuery(self, query_string, param_tuple)

    def show(self, table_name, max_field_size=20, html_table=False):
        """
        Show an entire table from the database.

        EXAMPLES:
            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_data('simon',[(0,0),(1,1),(1,2)])
            sage: DB.show('simon')
            a1                   b2
            ----------------------------------------
            0                    0
            1                    1
            1                    2
        """
        try:
            cur = self.__connection__.cursor()
            cur.execute('SELECT * FROM ' + table_name)
        except:
            raise RuntimeError('Failure to fetch data.')

        from sage.server.support import EMBEDDED_MODE
        if EMBEDDED_MODE or html_table:
            # Notebook Version
            print '<html><!--notruncate--><table bgcolor=lightgrey cellpadding=0><tr>'
            for des in cur.description:
                print '<td bgcolor=white align=center> ' + des[0] + ' </td>'
            print '</tr>'
            field_indices = range(len(cur.description))
            for row in cur:
                print '<tr>'
                for index in field_indices:
                    print '<td bgcolor=white align=center> ' + str(row[index]) + ' </td>'
                print '</tr>'
            print '</table></html>'

        else:
            # Command Prompt Version
            for des in cur.description:
                print des[0].ljust(max_field_size),
            print # new line
            print '-' * max_field_size * len(cur.description)
            field_indices = range(len(cur.description))
            for row in cur:
                for index in field_indices:
                    field_val = str(row[index])
                    print field_val.ljust(max_field_size) ,
                print # new line


class SQLDatabase(GenericSQLDatabase):

    def __init__(self, filename=None, skeleton=None):
        r"""
        A SQL Database object corresponding to a database file.

        INPUT:
            filename -- a string
            skeleton -- a triple-indexed dictionary
                    outer key - table name
                        inner key - column name
                            inner inner key - one of the following:
                    primary_key - boolean, whether column has been set as primary key
                    index - boolean, whether column has been set as index
                    sql - one of 'TEXT', 'BOOLEAN', 'INTEGER', 'REAL', or other
                        user defined type

        TUTORIAL:
        The SQLDatabase class is for interactively building databases intended for
        queries. This may sound redundant, but it is important. If you want a
        database intended for quick lookup of entries in very large tables, much
        like a hash table (such as a Python dictionary), a SQLDatabase may not be
        what you are looking for. The strength of SQLDatabases is in queries,
        searches through the database with complicated criteria.

        The class GenericSQLDatabase is for developers to provide a static
        database. The class does not support modification, and is meant to be a
        base class for specific classes of databases, such as the graph database.

        For example, we create a new database for storing isomorphism classes of
        simple graphs:
            sage: D = SQLDatabase()

        In order to generate representatives for the classes, we will import a
        function which generates all labeled graphs (noting that this is not the
        optimal way):
            sage: from sage.graphs.graph_isom import all_labeled_graphs

        We will need a table in the database in which to store the graphs, and we
        specify its structure with a Python dictionary, each of whose keys is the
        name of a column:
            sage: table_skeleton = {
            ... 'graph6':{'sql':'TEXT', 'index':True, 'primary_key':True},
            ... 'vertices':{'sql':'INTEGER'},
            ... 'edges':{'sql':'INTEGER'}
            ... }

        Then we create the table:
            sage: D.create_table('simon', table_skeleton)
            sage: D.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------

        Now that we have the table, we will begin to populate the table with
        rows. First, add the graph on zero vertices.
            sage: G = Graph()

            sage: D.add_row('simon',(0, G.graph6_string(), 0))

            sage: D.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------
            0                    ?                    0

        Next, add the graph on one vertex.
            sage: G.add_vertex()
            sage: D.add_row('simon',(0, G.graph6_string(), 1))
            sage: D.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------
            0                    ?                    0
            0                    @                    1

        Say we want a database of graphs on four or less vertices:
            sage: labels = {}
            sage: for i in range(2, 5):
            ...       labels[i] = []
            ...       for g in all_labeled_graphs(i):
            ...           g = g.canonical_label()
            ...           if g not in labels[i]:
            ...               labels[i].append(g)
            ...               D.add_row('simon', (g.size(), g.graph6_string(), g.order()))
            ...
            sage: D.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------
            0                    ?                    0
            0                    @                    1
            0                    A?                   2
            1                    A_                   2
            0                    B?                   3
            1                    BG                   3
            2                    BW                   3
            3                    Bw                   3
            0                    C?                   4
            1                    C@                   4
            2                    CB                   4
            3                    CF                   4
            3                    CJ                   4
            2                    CK                   4
            3                    CL                   4
            4                    CN                   4
            4                    C]                   4
            5                    C^                   4
            6                    C~                   4

        We can then query the database-- let's ask for all the graphs on four
        vertices with three edges. We do so by creating two queries and asking for
        rows that satisfy them both:
            sage: Q = SQLQuery(D, {'table_name':'simon', 'display_cols':['graph6'], 'expression':['vertices','=',4]})
            sage: Q2 = SQLQuery(D, {'table_name':'simon', 'display_cols':['graph6'], 'expression':['edges','=',3]})
            sage: Q = Q.intersect(Q2)
            sage: Q.run_query()
            [(u'CF', u'CF'), (u'CJ', u'CJ'), (u'CL', u'CL')]

        NOTE - The values of display_cols are always concatenated in intersections
        and unions.

        Of course, we can save the database to file:
            sage: replace_with_your_own_filepath = tmp_dir()
            sage: D.save(replace_with_your_own_filepath + 'simon.db')

        Now the database's hard link is to this file, and not the temporary db
        file. For example, let's say we open the same file with another class
        instance. We can load the file as an immutable database:
            sage: E = GenericSQLDatabase(replace_with_your_own_filepath + 'simon.db')
            sage: E.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------
            0                    ?                    0
            0                    @                    1
            0                    A?                   2
            1                    A_                   2
            0                    B?                   3
            1                    BG                   3
            2                    BW                   3
            3                    Bw                   3
            0                    C?                   4
            1                    C@                   4
            2                    CB                   4
            3                    CF                   4
            3                    CJ                   4
            2                    CK                   4
            3                    CL                   4
            4                    CN                   4
            4                    C]                   4
            5                    C^                   4
            6                    C~                   4
            sage: E.drop_table('simon')
            Traceback (most recent call last):
            ...
            AttributeError: 'GenericSQLDatabase' object has no attribute 'drop_table'
        """
        if filename is None:
            filename = tmp_filename() + '.db'
        elif (filename[-3:] != '.db'):
            raise ValueError('Please enter a valid database path (file name %s does not end in .db).'%filename)
        self.__dblocation__ = filename
        self.__connection__ = sqlite.connect(self.__dblocation__, check_same_thread=False)
        # check_same_thread = False:
        # this is to avoid the multiple thread problem with dsage:
        # pysqlite does not trust multiple threads for the same connection
        self.__connection__.create_function("regexp", 2, regexp)

        # construct skeleton (from provided database)
        self.__skeleton__ = construct_skeleton(self.__connection__)

        # add bones from new skeleton to database,
        # without changing existing structure
        if skeleton is not None:
            for table in skeleton:
                if table not in self.__skeleton__:
                    self.create_table(table, skeleton[table])
                else:
                    for column in skeleton[table]:
                        if column not in self.__skeleton__[table]:
                            self.create_column(table, column, skeleton[table][column])
                        else:
                            print "Column attributes were ignored for table %s, column %s -- column is already in table."%(table, column)

    def get_cursor(self):
        """
        Returns a pysqlite cursor for the database connection.

        A cursor is an input from which you can execute sqlite commands on the
        database.

        Recommended for more advanced users only.

        EXAMPLE:
            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_row('simon',(0,1))
            sage: DB.add_rows('simon',[(0,0),(1,1),(1,2)])
            sage: DB.add_rows('simon',[(0,0),(4,0),(5,1)], ['b2','a1'])
            sage: cur = DB.get_cursor()
            sage: (cur.execute('select * from simon')).fetchall()
            [(0, 1), (0, 0), (1, 1), (1, 2), (0, 0), (0, 4), (1, 5)]
        """
        return self.__connection__.cursor()

    def get_connection(self):
        """
        Returns a pysqlite connection to the database.

        You most likely want get_cursor() instead, which is used for executing
        sqlite commands on the database.

        Recommended for more advanced users only.

        EXAMPLE:
            sage: D = SQLDatabase()
            sage: con = D.get_connection()
            sage: t = con.execute('create table simon(n INTEGER, n2 INTEGER)')
            sage: for n in range(10):
            ...     t = con.execute('insert into simon values(%d,%d)'%(n,n^2))
            sage: D.show('simon')
            n                    n2
            ----------------------------------------
            0                    0
            1                    1
            2                    4
            3                    9
            4                    16
            5                    25
            6                    36
            7                    49
            8                    64
            9                    81
        """
        return self.__connection__

    def create_table(self, table_name, table_skeleton):
        """
        Creates a new table in the database.

        To create a table, a column structure must be specified. The form for
        this is a Python dict, for example:
        {'col1': {'sql':'INTEGER', 'index':False, 'primary_key':False}, ...}

        INPUT:
            table_name -- a string
            table_skeleton -- a double-indexed dictionary
                outer key - column name
                    inner key - one of the following:
                primary_key - boolean, whether column has been set as primary key
                index - boolean, whether column has been set as index
                sql - one of 'TEXT', 'BOOLEAN', 'INTEGER', 'REAL', or other
                    user defined type

        EXAMPLE:
            sage: D = SQLDatabase()
            sage: table_skeleton = {
            ... 'graph6':{'sql':'TEXT', 'index':True, 'primary_key':True},
            ... 'vertices':{'sql':'INTEGER'},
            ... 'edges':{'sql':'INTEGER'}
            ... }
            sage: D.create_table('simon', table_skeleton)
            sage: D.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------

        NOTE:
            Some SQL features, such as automatically incrementing primary key,
            require the full word 'INTEGER', not just 'INT'.

        """
        if self.__skeleton__.has_key(table_name):
            raise ValueError("Database already has a table named %s."%table_name)

        create_statement = 'create table ' + table_name + '('
        for col in table_skeleton:
            if col.find('sqlite') != -1:
                raise ValueError("Column names cannot contain 'sqlite'.")
            table_skeleton[col] = verify_column(table_skeleton[col])
            type = table_skeleton[col]['sql']
            if verify_type(type):
                if table_skeleton[col].has_key('primary_key') and table_skeleton[col]['primary_key']:
                    create_statement += col + ' ' + type + ' primary key, '
                else:
                    create_statement += col + ' ' + type + ', '
        create_statement = create_statement.rstrip(', ') + ') '

        self.__connection__.execute(create_statement)
        new_table_set_col_attr(self.__connection__, table_name, table_skeleton)
        self.__skeleton__[table_name] = table_skeleton

    def add_column(self, table_name, col_name, col_dict, default='NULL'):
        """
        Add a column named col_name to table table_name, whose data types are
        described by col_dict. The format for this is:
        {'col1':{'primary_key':False, 'index':True, 'sql':'REAL'}}

        INPUT:
            col_dict - a dictionary:
                key - column name
                    inner key - one of the following:
                primary_key - boolean, whether column has been set as primary key
                index - boolean, whether column has been set as index
                sql - one of 'TEXT', 'BOOLEAN', 'INTEGER', 'REAL', or other
                    user defined type

        EXAMPLES:
            sage: MonicPolys = SQLDatabase()
            sage: MonicPolys.create_table('simon', {'n':{'sql':'INTEGER', 'index':True}})
            sage: MonicPolys.show('simon')
            n
            --------------------
            sage: for n in range(20):
            ...       MonicPolys.add_row('simon', (n,))
            ...
            sage: MonicPolys.add_column('simon', 'n_squared', {'sql':'INTEGER', 'index':False}, 0)
            sage: MonicPolys.show('simon')
            n                    n_squared
            ----------------------------------------
            0                    0
            1                    0
            2                    0
            3                    0
            4                    0
            5                    0
            6                    0
            7                    0
            8                    0
            9                    0
            10                   0
            11                   0
            12                   0
            13                   0
            14                   0
            15                   0
            16                   0
            17                   0
            18                   0
            19                   0
            sage: MonicPolys.drop_column('simon', 'n_squared')
            sage: MonicPolys.show('simon')
            n
            --------------------
            0
            1
            2
            3
            4
            5
            6
            7
            8
            9
            10
            11
            12
            13
            14
            15
            16
            17
            18
            19

        """
        # Check input:
        if col_name.find('sqlite') != -1:
            raise ValueError("Column names cannot contain 'sqlite'.")
        if not self.__skeleton__.has_key(table_name):
            raise ValueError("Database has no table %s."%table_name)
        if self.__skeleton__[table_name].has_key(col_name):
            raise ValueError("Table %s already has column %s."%(table_name,col_name))
        col_dict = verify_column(col_dict)

        # Get an ordered list:
        cur_list = skel_to_col_attr_list(self.__skeleton__[table_name])
        # Update the skeleton:
        self.__skeleton__[table_name][col_name] = col_dict

        original = ''
        for col in cur_list:
            original += col[0] +', '
        original = original.rstrip(', ')

        more = original + ', ' + col_name
        more_attr = ''
        for col in cur_list:
            if col[2]: # If primary key:
                more_attr += col[0] + ' ' + col[1] + ' primary key, '
            else:
                more_attr += col[0] + ' ' + col[1] + ', '
        more_attr += col_name + ' ' + col_dict['sql']
        try:
            # Silly SQLite -- we have to make a temp table to hold info...
            self.__connection__.executescript("""
                create temporary table spam(%s);
                insert into spam select %s, %s from %s;
                drop table %s;
                create table %s (%s);
                """%(more_attr, original, default, table_name, table_name, table_name, more_attr))

            # Update indices in new table
            new_table_set_col_attr(self.__connection__, table_name, self.__skeleton__[table_name])

            # Now we can plop our data into the *new* table:
            self.__connection__.executescript("""
                insert into %s select %s from spam;
                drop table spam;
                """%(table_name, more))

            self.vacuum()
        except sqlite.Error, e:
            print 'A sqlite error occured: ', e.args[0]
            # delete added column from skeleton
            self.__skeleton__[table_name].pop(col_name)

    def drop_column(self, table_name, col_name):
        """
        Drop the column col_name from table table_name.

        EXAMPLES:
            sage: MonicPolys = SQLDatabase()
            sage: MonicPolys.create_table('simon', {'n':{'sql':'INTEGER', 'index':True}})
            sage: MonicPolys.show('simon')
            n
            --------------------
            sage: for n in range(20):
            ...       MonicPolys.add_row('simon', (n,))
            ...
            sage: MonicPolys.add_column('simon', 'n_squared', {'sql':'INTEGER', 'index':False}, 0)
            sage: MonicPolys.show('simon')
            n                    n_squared
            ----------------------------------------
            0                    0
            1                    0
            2                    0
            3                    0
            4                    0
            5                    0
            6                    0
            7                    0
            8                    0
            9                    0
            10                   0
            11                   0
            12                   0
            13                   0
            14                   0
            15                   0
            16                   0
            17                   0
            18                   0
            19                   0
            sage: MonicPolys.drop_column('simon', 'n_squared')
            sage: MonicPolys.show('simon')
            n
            --------------------
            0
            1
            2
            3
            4
            5
            6
            7
            8
            9
            10
            11
            12
            13
            14
            15
            16
            17
            18
            19

        """
        # Check input:
        if not self.__skeleton__.has_key(table_name):
            raise ValueError("Database has no table %s."%table_name)
        if not self.__skeleton__[table_name].has_key(col_name):
            raise ValueError("Table %s has no column %s."%(table_name,col_name))

        # Update the skeleton:
        self.__skeleton__[table_name].pop(col_name)
        # Get an ordered list (without the column we're deleting):
        cur_list = skel_to_col_attr_list(self.__skeleton__[table_name])

        less = ''
        for col in cur_list:
            less += col[0] +', '
        less = less.rstrip(', ')

        less_attr = ''
        less_attr = ''
        for col in cur_list:
            if col[2]: # If primary key:
                less_attr += col[0] + ' ' + col[1] + ' primary key, '
            else:
                less_attr += col[0] + ' ' + col[1] + ', '
        less_attr = less_attr.rstrip(', ')

        # Silly SQLite -- we have to make a temp table to hold info...
        self.__connection__.executescript("""
            create temporary table spam(%s);
            insert into spam select %s from %s;
            drop table %s;
            create table %s (%s);
            """%(less_attr, less, table_name, table_name, table_name, less_attr))
        # Update indices in new table
        new_table_set_col_attr(self.__connection__, table_name, self.__skeleton__[table_name])

        # Now we can plop our data into the *new* table:
        self.__connection__.executescript("""
            insert into %s select %s from spam;
            drop table spam;
            """%(table_name, less))

        self.vacuum()

    def rename_table(self, table_name, new_name):
        """
        Renames the table table_name to new_name.

        EXAMPLE:
            sage: D = SQLDatabase()
            sage: D.create_table('simon',{'col1':{'sql':'INTEGER'}})
            sage: D.show('simon')
            col1
            --------------------
            sage: D.rename_table('simon', 'lucy')
            sage: D.show('simon')
            Traceback (most recent call last):
            ...
            RuntimeError: Failure to fetch data.

            sage: D.show('lucy')
            col1
            --------------------

        """
        # Check input:
        if not self.__skeleton__.has_key(table_name):
            raise ValueError("Database has no table %s."%table_name)
        if self.__skeleton__.has_key(new_name):
            raise ValueError("Database already has table %s."%new_name)

        self.__connection__.execute('alter table %s rename to %s'%(table_name, new_name))

        # Update skeleton:
        self.__skeleton__[new_name] = self.__skeleton__[table_name]
        self.__skeleton__.pop(table_name)

    def drop_table(self, table_name):
        """
        Delete table table_name from database.

        INPUT:
            table_name -- a string

        EXAMPLE:
            sage: D = SQLDatabase()
            sage: D.create_table('simon',{'col1':{'sql':'INTEGER'}})
            sage: D.show('simon')
            col1
            --------------------
            sage: D.drop_table('simon')
            sage: D.get_skeleton()
            {}

        """
        if not self.__skeleton__.has_key(table_name):
            raise ValueError("Database has no table %s."%table_name)

        self.__connection__.execute('drop table ' + table_name)

        # Update Skeleton
        self.__skeleton__.pop(table_name)

    def drop_data_from_table(self, table_name):
        """
        Removes all data from table_name, except for the structure of the
        columns.

        EXAMPLE:
            sage: D = SQLDatabase()
            sage: D.create_table('simon',{'col1':{'sql':'INTEGER'}})
            sage: D.add_row('simon',(9,))
            sage: D.show('simon')
            col1
            --------------------
            9
            sage: D.drop_data_from_table('simon')
            sage: D.show('simon')
            col1
            --------------------

        """
        if not self.__skeleton__.has_key(table_name):
            raise ValueError("Database has no table %s."%table_name)
        self.__connection__.execute('delete from ' + table_name)

    def make_index(self, col_name, table_name, unique=False):
        """
        Set the column col_name in table table_name to be an index, that is, a
        column set up to do quick searches on.

        INPUT:
            col_name -- a string
            table_name -- a string
            unique -- requires that there are no multiple entries in the
                column, makes searching faster

        EXAMPLE:
            sage: MonicPolys = SQLDatabase()
            sage: MonicPolys.create_table('simon', {'n':{'sql':'INTEGER', 'index':True}, 'n2':{'sql':'INTEGER'}})
            sage: for n in range(20):
            ...       MonicPolys.add_row('simon', (n**2,n))
            ...
            sage: MonicPolys.show('simon')
            n2                   n
            ----------------------------------------
            0                    0
            1                    1
            4                    2
            9                    3
            16                   4
            25                   5
            36                   6
            49                   7
            64                   8
            81                   9
            100                  10
            121                  11
            144                  12
            169                  13
            196                  14
            225                  15
            256                  16
            289                  17
            324                  18
            361                  19
            sage: MonicPolys.make_index('n2','simon')
            sage: MonicPolys.get_skeleton()
            {'simon': {'n2': {'index': True, 'primary_key': False, 'sql': 'INTEGER'}, 'n': {'index': True, 'primary_key': False, 'sql': 'INTEGER'}}}
            sage: MonicPolys.drop_index('simon', 'n')
            sage: MonicPolys.make_primary_key('simon', 'n2')
            sage: MonicPolys.get_skeleton()
            {'simon': {'n2': {'index': True, 'primary_key': True, 'sql': 'INTEGER'}, 'n': {'index': False, 'primary_key': False, 'sql': 'INTEGER'}}}
            sage: MonicPolys.drop_primary_key('simon', 'n2')
            sage: MonicPolys.get_skeleton()
            {'simon': {'n2': {'index': True, 'primary_key': False, 'sql': 'INTEGER'}, 'n': {'index': False, 'primary_key': False, 'sql': 'INTEGER'}}}
        """
        if not self.__skeleton__.has_key(table_name):
            raise ValueError("Database has no table %s."%table_name)
        if not self.__skeleton__[table_name].has_key(col_name):
            raise ValueError("Table %s has no column %s."%(table_name,col_name))

        if unique:
            index_string = 'create unique index ' + col_name + ' on ' + table_name + ' (' + col_name + ')'
        else:
            index_string = 'create index ' + col_name + ' on ' + table_name + ' (' + col_name + ')'
        cur = self.__connection__.cursor()
        exe = cur.execute(index_string)

        # Update Skeleton
        self.__skeleton__[table_name][col_name]['index'] = True

    def drop_index(self, table_name, index_name):
        """
        Set the column index_name in table table_name to not be an index. See
        make_index()

        EXAMPLE:
            sage: MonicPolys = SQLDatabase()
            sage: MonicPolys.create_table('simon', {'n':{'sql':'INTEGER', 'index':True}, 'n2':{'sql':'INTEGER'}})
            sage: for n in range(20):
            ...       MonicPolys.add_row('simon', (n**2,n))
            ...
            sage: MonicPolys.show('simon')
            n2                   n
            ----------------------------------------
            0                    0
            1                    1
            4                    2
            9                    3
            16                   4
            25                   5
            36                   6
            49                   7
            64                   8
            81                   9
            100                  10
            121                  11
            144                  12
            169                  13
            196                  14
            225                  15
            256                  16
            289                  17
            324                  18
            361                  19
            sage: MonicPolys.make_index('n2','simon')
            sage: MonicPolys.get_skeleton()
            {'simon': {'n2': {'index': True, 'primary_key': False, 'sql': 'INTEGER'}, 'n': {'index': True, 'primary_key': False, 'sql': 'INTEGER'}}}
            sage: MonicPolys.drop_index('simon', 'n')
            sage: MonicPolys.make_primary_key('simon', 'n2')
            sage: MonicPolys.get_skeleton()
            {'simon': {'n2': {'index': True, 'primary_key': True, 'sql': 'INTEGER'}, 'n': {'index': False, 'primary_key': False, 'sql': 'INTEGER'}}}
            sage: MonicPolys.drop_primary_key('simon', 'n2')
            sage: MonicPolys.get_skeleton()
            {'simon': {'n2': {'index': True, 'primary_key': False, 'sql': 'INTEGER'}, 'n': {'index': False, 'primary_key': False, 'sql': 'INTEGER'}}}

        """
        if not self.__skeleton__.has_key(table_name):
            raise ValueError("Database has no table %s."%table_name)
        if not self.__skeleton__[table_name].has_key(index_name):
            raise ValueError("Table %s has no column %s."%(table,index_name))
        if not self.__skeleton__[table_name][index_name]['index']:
            return # silently

        cur = self.__connection__.cursor()
        exe = cur.execute('drop index ' + index_name)

        # Update Skeleton
        self.__skeleton__[table_name][index_name]['index'] = False

    def make_primary_key(self, table_name, col_name):
        """
        Set the column col_name in table table_name to be a primary key.

        A primary key is something like an index, but its main purpose is to
        link different tables together. This allows searches to be executed on
        multiple tables that represent maybe different data about the same
        objects.

        NOTE:
            Some SQL features, such as automatically incrementing primary key,
            require the full word 'INTEGER', not just 'INT'.

        EXAMPLE:
            sage: MonicPolys = SQLDatabase()
            sage: MonicPolys.create_table('simon', {'n':{'sql':'INTEGER', 'index':True}, 'n2':{'sql':'INTEGER'}})
            sage: for n in range(20):
            ...       MonicPolys.add_row('simon', (n**2,n))
            ...
            sage: MonicPolys.show('simon')
            n2                   n
            ----------------------------------------
            0                    0
            1                    1
            4                    2
            9                    3
            16                   4
            25                   5
            36                   6
            49                   7
            64                   8
            81                   9
            100                  10
            121                  11
            144                  12
            169                  13
            196                  14
            225                  15
            256                  16
            289                  17
            324                  18
            361                  19
            sage: MonicPolys.make_index('n2','simon')
            sage: MonicPolys.get_skeleton()
            {'simon': {'n2': {'index': True, 'primary_key': False, 'sql': 'INTEGER'}, 'n': {'index': True, 'primary_key': False, 'sql': 'INTEGER'}}}
            sage: MonicPolys.drop_index('simon', 'n')
            sage: MonicPolys.make_primary_key('simon', 'n2')
            sage: MonicPolys.get_skeleton()
            {'simon': {'n2': {'index': True, 'primary_key': True, 'sql': 'INTEGER'}, 'n': {'index': False, 'primary_key': False, 'sql': 'INTEGER'}}}
            sage: MonicPolys.drop_primary_key('simon', 'n2')
            sage: MonicPolys.get_skeleton()
            {'simon': {'n2': {'index': True, 'primary_key': False, 'sql': 'INTEGER'}, 'n': {'index': False, 'primary_key': False, 'sql': 'INTEGER'}}}
        """

        # To get around sqlite's problems with *alter table* commands, we
        # have set up the creation of a temporary database in this method.
        # Please feel free to improve on this by sending a patch or suggestion.

        if not self.__skeleton__.has_key(table_name):
            raise ValueError("Database has no table %s."%table_name)
        if not self.__skeleton__[table_name].has_key(col_name):
            raise ValueError("Table %s has no column %s."%(table_name,col_name))

        # Update the skeleton:
        self.__skeleton__[table_name][col_name]['primary_key'] = True
        # Get an ordered list (with the primary key info updated):
        cur_list = skel_to_col_attr_list(self.__skeleton__[table_name])

        new = ''
        for col in cur_list:
            new += col[0] +', '
        new = new.rstrip(', ')

        new_attr = ''
        new_attr = ''
        for col in cur_list:
            if col[2]: # If primary key:
                new_attr += col[0] + ' ' + col[1] + ' primary key, '
            else:
                new_attr += col[0] + ' ' + col[1] + ', '
        new_attr = new_attr.rstrip(', ')

        # Silly SQLite -- we have to make a temp table to hold info...
        self.__connection__.executescript("""
            create temporary table spam(%s);
            insert into spam select %s from %s;
            drop table %s;
            create table %s (%s);
            """%(new_attr, new, table_name, table_name, table_name,new_attr))

        # Update indices in new table
        new_table_set_col_attr(self.__connection__, table_name, self.__skeleton__[table_name])

        # Now we can plop our data into the *new* table:
        self.__connection__.executescript("""
            insert into %s select %s from spam;
            drop table spam;
            """%(table_name, new))
        self.vacuum()

    def drop_primary_key(self, table_name, col_name):
        """
        Set the column col_name in table table_name not to be a primary key.

        A primary key is something like an index, but its main purpose is to
        link different tables together. This allows searches to be executed on
        multiple tables that represent maybe different data about the same
        objects.

        Note: This function only changes the column to be non-primary, it does
        not delete it.

        EXAMPLE:
            sage: MonicPolys = SQLDatabase()
            sage: MonicPolys.create_table('simon', {'n':{'sql':'INTEGER', 'index':True}, 'n2':{'sql':'INTEGER'}})
            sage: for n in range(20):
            ...       MonicPolys.add_row('simon', (n**2,n))
            ...
            sage: MonicPolys.show('simon')
            n2                   n
            ----------------------------------------
            0                    0
            1                    1
            4                    2
            9                    3
            16                   4
            25                   5
            36                   6
            49                   7
            64                   8
            81                   9
            100                  10
            121                  11
            144                  12
            169                  13
            196                  14
            225                  15
            256                  16
            289                  17
            324                  18
            361                  19
            sage: MonicPolys.make_index('n2','simon')
            sage: MonicPolys.get_skeleton()
            {'simon': {'n2': {'index': True, 'primary_key': False, 'sql': 'INTEGER'}, 'n': {'index': True, 'primary_key': False, 'sql': 'INTEGER'}}}
            sage: MonicPolys.drop_index('simon', 'n')
            sage: MonicPolys.make_primary_key('simon', 'n2')
            sage: MonicPolys.get_skeleton()
            {'simon': {'n2': {'index': True, 'primary_key': True, 'sql': 'INTEGER'}, 'n': {'index': False, 'primary_key': False, 'sql': 'INTEGER'}}}
            sage: MonicPolys.drop_primary_key('simon', 'n2')
            sage: MonicPolys.get_skeleton()
            {'simon': {'n2': {'index': True, 'primary_key': False, 'sql': 'INTEGER'}, 'n': {'index': False, 'primary_key': False, 'sql': 'INTEGER'}}}

        """
        # To get around sqlite's problems with *alter table* commands, we
        # have set up the creation of a temporary database in this method.
        # Please feel free to improve on this by sending a patch or suggestion.

        if not self.__skeleton__.has_key(table_name):
            raise ValueError("Database has no table %s."%table_name)
        if not self.__skeleton__[table_name].has_key(col_name):
            raise ValueError("Table %s has no column %s."%(table_name,col_name))
        if not self.__skeleton__[table_name][col_name]['primary_key']:
            return # silently

        # Update the skeleton:
        self.__skeleton__[table_name][col_name]['primary_key'] = False
        # Get an ordered list (with the primary key info updated):
        cur_list = skel_to_col_attr_list(self.__skeleton__[table_name])

        new = ''
        for col in cur_list:
            new += col[0] +', '
        new = new.rstrip(', ')

        new_attr = ''
        new_attr = ''
        for col in cur_list:
            if col[2]: # If primary key:
                new_attr += col[0] + ' ' + col[1] + ' primary key, '
            else:
                new_attr += col[0] + ' ' + col[1] + ', '
        new_attr = new_attr.rstrip(', ')

        # Silly SQLite -- we have to make a temp table to hold info...
        self.__connection__.executescript("""
            create temporary table spam(%s);
            insert into spam select %s from %s;
            drop table %s;
            create table %s (%s);
            """%(new_attr, new, table_name, table_name, table_name, new_attr))

        # Update indices in new table
        new_table_set_col_attr(self.__connection__, table_name, self.__skeleton__[table_name])

        # Now we can plop our data into the *new* table:
        self.__connection__.executescript("""
            insert into %s select %s from spam;
            drop table spam;
            """%(table_name, new))
        self.vacuum()

    def add_row(self, table_name, values, entry_order=None):
        """
        Add the row described by values to the table table_name. Values should
        be a tuple, of same length and order as columns in given table.

        NOTE: If values is of length one, be sure to specify that it is a tuple
        of length one, by using a comma, e.g.:
            sage: values = (6,)

        EXAMPLES:
            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_row('simon',(0,1))
            sage: cur = DB.get_cursor()
            sage: (cur.execute('select * from simon')).fetchall()
            [(0, 1)]
        """
        self.add_data(table_name, [values], entry_order)

    def delete_rows(self, query):
        """
        Uses a SQLQuery instance to modify (delete rows from) the
        database.  Note that this function will not allow deletion via a
        GenericSQLQuery (a method for more advanced users) in order to
        prevent an accidental disaster (omitting a where clause or using '*').

        SQLQuery must have no join statements.  (As of now, you can only
        delete from one table at a time -- ideas and patches are welcome).

        To remove all data that satisfies a SQLQuery, send the query as an
        argument to delete_rows.  Be careful, test your query first.

        Recommended use:  have some kind of row identification primary
        key column that you use as a parameter in the query.  (See example
        below).

        INPUT:
            query -- a SQLQuery (Delete the rows returned when query is
                run).

        EXAMPLES:
            sage: DB = SQLDatabase()
            sage: DB.create_table('lucy',{'id':{'sql':'INTEGER', 'primary_key':True, 'index':True},'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_rows('lucy', [(0,1,1),(1,1,4),(2,0,7),(3,1,384),(4,1,978932)],['id','a1','b2'])
            sage: DB.show('lucy')
            a1                   id                   b2
            ------------------------------------------------------------
            1                    0                    1
            1                    1                    4
            0                    2                    7
            1                    3                    384
            1                    4                    978932
            sage: Q = SQLQuery(DB, {'table_name':'lucy', 'display_cols':['id','a1','b2'], 'expression':['id','>=',3]})
            sage: Q.show()
            id                   a1                   b2
            ------------------------------------------------------------
            3                    1                    384
            4                    1                    978932
            sage: DB.delete_rows(Q)
            sage: DB.show('lucy')
            a1                   id                   b2
            ------------------------------------------------------------
            1                    0                    1
            1                    1                    4
            0                    2                    7

        """
        # Check query is associated with this database
        if not isinstance(query, SQLQuery):
            raise TypeError('%s is not a valid SQLQuery'%query)
        if query.__database__ is not self:
            raise ValueError('%s is not associated to this database.'%query)
        if (query.__query_string__).find(' JOIN ') != -1:
            raise ValueError('%s is not a valid query.  Can only delete from one table at a time.'%query)

        delete_statement = re.sub('SELECT .* FROM', 'DELETE FROM', query.__query_string__)

        try:
            tup = str(query.__param_tuple__).rstrip(')') + ',)'
            cur = self.__connection__.cursor()
            if query.__param_tuple__ is not None:
                tup = []
                for i in range(len(query.__param_tuple__)):
                    tup.append(str(query.__param_tuple__[i]))
                cur.execute(delete_statement, tuple(tup))
            else:
                cur.execute(delete_statement)
        except:
            raise RuntimeError('Failure to complete delete. Check your data.')

    def add_rows(self, table_name, rows, entry_order=None):
        """
        INPUT:
            rows is a list of tuples that represent one row of data to add
            (types should match col types in order)
            entry_order --  an ordered list or tuple
                            overrides normal order with user defined order

        EXAMPLES:
            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_row('simon',(0,1))
            sage: DB.add_rows('simon',[(0,0),(1,1),(1,2)])
            sage: DB.add_rows('simon',[(0,0),(4,0),(5,1)], ['b2','a1'])
            sage: cur = DB.get_cursor()
            sage: (cur.execute('select * from simon')).fetchall()
            [(0, 1), (0, 0), (1, 1), (1, 2), (0, 0), (0, 4), (1, 5)]
        """
        self.add_data(table_name, rows, entry_order)

    def add_data(self, table_name, rows, entry_order=None):
        """
        Add data from a list of rows to the database.

        INPUT:
            rows is a list of tuples that represent one row of data to add
            (types should match col types in order)
            entry_order --  an ordered list or tuple
                            overrides normal order with user defined order

        EXAMPLES:
            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_row('simon',(0,1))
            sage: DB.add_data('simon',[(0,0),(1,1),(1,2)])
            sage: DB.add_data('simon',[(0,0),(4,0),(5,1)], ['b2','a1'])
            sage: cur = DB.get_cursor()
            sage: (cur.execute('select * from simon')).fetchall()
            [(0, 1), (0, 0), (1, 1), (1, 2), (0, 0), (0, 4), (1, 5)]
        """
        quest = '('
        length = len(rows[0])
        for i in range (length):
            quest += '?, '
        quest = quest.rstrip(', ') + ')'
        strows = []
        for row in rows:
            tup = []
            for entry in row:
                tup.append(str(entry))
            strows.append(tuple(tup))

        if entry_order is not None:
            self.__connection__.executemany("INSERT INTO " + table_name + str(tuple(entry_order)) + " VALUES " + quest, strows)
        else:
            self.__connection__.executemany("INSERT INTO " + table_name + " VALUES " + quest, strows)

    def vacuum(self):
        """
        Cleans the extra hard disk space used up by a database that has
        recently shrunk.

        EXAMPLE:
            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_row('simon',(0,1))
            sage: DB.add_data('simon',[(0,0),(1,1),(1,2)])
            sage: DB.add_data('simon',[(0,0),(4,0),(5,1)], ['b2','a1'])
            sage: DB.show('simon')
            a1                   b2
            ----------------------------------------
            0                    1
            0                    0
            1                    1
            1                    2
            0                    0
            0                    4
            1                    5
            sage: DB.drop_column('simon','b2')
            sage: DB.commit()
            sage: DB.vacuum()
            sage: DB.show('simon')
            a1
            --------------------
            0
            0
            1
            1
            0
            0
            1
        """
        self.__connection__.execute('vacuum')

    def commit(self):
        """
        Commits changes to file.

        EXAMPLE:
            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_row('simon',(0,1))
            sage: DB.add_data('simon',[(0,0),(1,1),(1,2)])
            sage: DB.add_data('simon',[(0,0),(4,0),(5,1)], ['b2','a1'])
            sage: DB.show('simon')
            a1                   b2
            ----------------------------------------
            0                    1
            0                    0
            1                    1
            1                    2
            0                    0
            0                    4
            1                    5
            sage: DB.drop_column('simon','b2')
            sage: DB.commit()
            sage: DB.vacuum()
            sage: DB.show('simon')
            a1
            --------------------
            0
            0
            1
            1
            0
            0
            1
        """
        try:
            self.__connection__.execute('commit')
        except:
            # Not sure why this throws an exception - but without it,
            #       the changes are not committed so it is necessary.
            pass
