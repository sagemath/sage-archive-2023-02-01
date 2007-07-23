"""
TODO
 - docstrings and e.g.'s - (there is an example in query intersection to get you going)
 - There might be some nonetype errors throughout that should be easily
   found/fixed when making doctests.
 - Also, some query_dict (think, skeleton) checking on intersect and
   union hasn't been done yet.

EMILY: This is a new list:
1 - There should probably be an add_rows function, in the plural, where values
is a list of tuples...
2 - There is a bug at the bottom of the file I'm not exactly sure how to
handle. I think a try/except clause, but to do that we need to know what is
causing this behavior, and what effects it has, so we can undo them in the
except part.
3 - We need to figure out how temp files are handled, and if there is a way to
delete the database file as soon as the variable gets deallocated, if it is a
temp database file.
4 - If you look at the docstrings for SQLDatabase.add_column(), you might
notice that there is no way to now set the values in the new column, without
deleting all the rows. Maybe this is a SQLite thing...
5 - Documentation for index and primary key stuff.
6 - As far as try/except things, I only think we should add them in at this
point if there is something bad SQL is doing before it fails (see add_data).

Databases.


            skeleton -- a triple-indexed dictionary
                outer key - table name
                    inner key - column name
                        inner inner key - one of the following:
                primary_key - boolean, whether column has been set as primary key
                index - boolean, whether column has been set as index
                sql - one of 'STRING', 'BOOLEAN', 'INTEGER', 'REAL', or other
                    user defined type

        An example skeleton of a database with one table, that table with one
        column:
        {'table1':{'col1':{'primary_key':False, 'index':True, 'sql':'REAL'}}}

FUTURE TODOs (Ignore for now):
    - order by clause in query strings
    - delete from query containing joins
"""

################################################################################
#           Copyright (C) 2007 Emily A. Kirkman
#                              Robert L. Miller
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
################################################################################
from sqlite3 import dbapi2 as sqlite # TODO: maybe some comment why dbapi2 instead of sqlite3?
import os
import re
from sage.misc.misc import tmp_filename
from sage.structure.sage_object import SageObject
from sage.server.support import EMBEDDED_MODE

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

def verify_type(type):
    types = ['INTEGER','INT','BOOLEAN','REAL','STRING','BOOL']
    if type.upper() not in types:
        raise TypeError('%s is not a legal type.'%type)
    return True

def verify_column(col_dict):
    """
    Verify that a column dict is in proper format*, and return a dict with
    default values filled in.

    * {'primary_key':False, 'index':False, 'sql':'REAL'}

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
    binaries = ['=','<=','>=','like','<','>','<>','regexp']
    unaries = ['is null','is not null']
    if operator not in binaries and operator not in unaries:
        raise TypeError('%s is not a legal operator.'%operator)
    return True

def construct_skeleton(connection):
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
            skeleton[table[0]][column[1]]['index'] = True
    return skeleton

def skel_to_col_attr_list(table_dict):
    s = []
    for col in table_dict:
        s.append((col, table_dict[col]['sql'], table_dict[col]['primary_key']))
    return s

def new_table_set_col_attr(connection, table_name, table_skeleton):
    statement = ''
    for col in table_skeleton:
        if table_skeleton[col].has_key('index'):
            if table_skeleton[col]['index']:
                statement += 'CREATE INDEX %s ON %s (%s) '%(col, table_name, col)
        else:
            table_skeleton[col]['index'] = False
    connection.execute(statement)

class GenericSQLQuery(SageObject):
    """
    Emily - documentation

    ABSOLUTELY NO RESPONSIBILITY FOR THIS!!!
    """

    def __init__(self, database, query_string, param_tuple=None):
        """
        Emily - documentation

        TEACH PEOPLE ABOUT USING '?' AND TUPLES
        """

        if not isinstance(database, SQLDatabase):
            raise TypeError('%s is not a valid SQLDatabase'%database)

        self.__database__ = database
        self.__param_tuple__ = param_tuple
        self.__query_string__ = query_string

    def __repr__(self):
        """
        __repr__ gets called when you type self and hit enter. It should
        return a string representing the object. Here, the current query
        string along with the parameter tuples are printed.

        """
        s = "Query for sql database: "
        s += self.__database__ + "\n"
        s += "Query string: "
        s += self.__query_string__ + "\n"
        s += "Parameter tuple: "
        s += str(self.__param_tuple__) + "\n"
        return s

    def copy(self):
        """
        Emily - documentation

        """
        return GenericSQLQuery(self.__database__, self.__query_string__, self.__param_tuple__)

    def run_query(self):
        """
        Emily - documentation

        """
        try:
            #tup = str(self.__param_tuple__).rstrip(')') + ',)'
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

    def show(self, max_field_size=20):
        """
        Displays the result of the query in table format.

        INPUT:
            max_field_size -- how wide each field can be

        EXAMPLE:
            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_data('simon',[(0,0),(1,1),(1,2)])
            sage: r = SQLQuery(DB, {'table_name':'simon', 'display_cols':'a1', 'expression':['b2','<=', 6]})
            sage: p = SQLQuery(DB, {'table_name':'simon', 'display_cols':'b2', 'expression':['b2','<=', 6]})
            sage: s = p.intersect(r)
            sage: s.show()
            b2                   a1
            ----------------------------------------
            0                    0
            1                    1
            2                    1

        """
        try:
            #tup = str(self.__param_tuple__).rstrip(')') + ',)'
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

        if EMBEDDED_MODE:
            # Notebook Version
            print '<html><table bgcolor=lightgrey cellpadding=0><tr>'
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

class SQLQuery(GenericSQLQuery):
    """
    EMILY - documentation

    query_dict := {'table_name': 'tblname', 'display_cols': 'col1, col2, col3', 'expression':[col, operator, value]}

    point out strings '"value"'

    """

    def __init__(self, database, query_dict=None):
        if not isinstance(database, GenericSQLDatabase):
            raise TypeError('%s is not a valid SQLDatabase'%database)

        self.__database__ = database
        self.__query_dict__ = query_dict

        # TODO: ERROR CHECKING:
        # 1. Confirm query_dict matches skeleton:
        # 2. confirm tblname is in database:
        # 3. confirm display_cols in tblname:
        # 4. confirm col (from expression is in tblname):

        # confirm operator:
        verify_operator(query_dict['expression'][1])

        # make tuple:
        if self.__query_dict__ is not None:
            self.__param_tuple__ = (self.__query_dict__['expression'][2],)
        else:
            self.__param_tuple__ = None

        # make query string:
        if self.__query_dict__ is not None:
            self.__query_string__ = 'SELECT ' + self.__query_dict__['table_name'] + \
                                    '.' + self.__query_dict__['display_cols'] + \
                                    ' FROM ' + self.__query_dict__['table_name'] + \
                                    ' WHERE ' + self.__query_dict__['table_name'] + '.' + \
                                    self.__query_dict__['expression'][0] + ' ' + \
                                    self.__query_dict__['expression'][1] + ' ?'
        else:
            self.__query_string__ = None

    def copy(self):
        return SQLQuery(self.__database__, self.__query_dict__)

    def intersect(self, other, join_table=None, join_dict=None):
        """
        join_dict -- {join_table1: (corr_base_col1, col1), join_table2: (corr_base_col2, col2)}
        join_table -- base table to join on

        EXAMPLES:
            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.create_table('lucy',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_data('simon', [(0,5),(1,4)])
            sage: DB.add_data('lucy', [(1,1),(1,4)])
            sage: q = SQLQuery(DB, {'table_name':'lucy', 'display_cols':'b2', 'expression':['a1','=',1]})
            sage: r = SQLQuery(DB, {'table_name':'simon', 'display_cols':'a1', 'expression':['b2','<=', 6]})
            sage: s = q.intersect(r, 'simon', {'lucy':('a1','a1')})
            sage: s.__query_string__
            'SELECT lucy.b2,simon.a1 FROM simon INNER JOIN lucy ON simon.a1=lucy.a1 WHERE ( lucy.a1 = ? ) AND ( simon.b2 <= ? )'
            sage: s.run_query()
            [(1, 1), (4, 1)]
        """

        # TODO : Check same database for all joins (i.e.: join_table1 and join_table2 are in same database)

        # TODO : SOME CHECKING - (1) if more than one table (AT ALL INVOLVED WHATSOEVER), both join args must not be None
        #                        (2) also, check in with (database) skeleton for possible bad input
        #    NOT ROBERT-         (3) and compare old query strings to confirm all previous
        #                           tables are included - (otherwise Error)
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

        q.__param_tuple__ = self.__param_tuple__ + other.__param_tuple__

        return q

    def union(self, other, join_table=None, join_dict=None):
        """
        join_dict -- {join_table1: (corr_base_col1, col1), join_table2: (corr_base_col2, col2)}
        """

        # TODO : Check same database for all joins (i.e.: join_table1 and join_table2 are in same database)

        # TODO : SOME CHECKING - (1) if more than one table, both join args must not be None
        #                    (2) also, check in with skeleton for possible bad input
        #                    (3) and compare old query strings to confirm all previous
        #                           tables are included - (otherwise Error)
        q = self.copy()

        # inner join clause
        if join_dict is not None:
            joins = join_table
            for table in join_dict:
                joins += ' INNER JOIN ' + table + ' ON ' + join_table + \
                         '.' + join_dict[table][0] + '=' + table + '.' + \
                         join_dict[table][1] + ' '
            q.__query_string__ = re.sub(' FROM .* WHERE ', ' FROM' + joins + 'WHERE ', self.__query_string__)

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

    def complement(self):
        q = SQLQuery(self.__database__)
        q.__query_string__ = re.sub(' WHERE ',' WHERE NOT ( ',self.__query_string__)
        q.__query_string__ += ' )'
        q.__param_tuple__ = self.__param_tuple__
        return q

class GenericSQLDatabase(SageObject):
    """
    *Immutable* Database class.

    INPUT:
        filename -- where to load the database from

    TODO:
        once the tutorial is finished in SQLDatabase, copy it here

    """
    def __init__(self, filename):

        if (filename[-3:] != '.db'):
            raise ValueError('Please enter a valid database path (file name %s does not end in .db).'%filename)
        self.__dblocation__ = filename
        self.__connection__ = sqlite.connect(self.__dblocation__)

        self.__skeleton__ = construct_skeleton(self.__connection__)

    def __repr__(self):
        s = ''
        for table in self.__skeleton__:
            s += 'table ' + table + ':\n'
            for column in self.__skeleton__[table]:
                s += '   column ' + column + ': '
                for data in self.__skeleton__[table][column]:
                    s += data + ': ' + self.__skeleton__[table][column][data] + '; '
                s += '\n'
        return s

    def copy(self):
        """
        Returns an instance of SQLDatabase that points to a copy database,
        and allows modification.

        TODO - examples (once there are immutable databases to start from)
        """
        # copy .db file
        new_loc = tmp_filename() + '.db'
        os.system('cp '+ self.__dblocation__ + ' ' + new_loc)
        D = SQLDatabase(filename=new_loc, skeleton=copy(self.__skeleton__))
        return D

    def save(self, filename):
        """
        Save the database to the specified location.

        TODO - figure out bug regarding saving then loading, and use something
        similar to the example in the tutorial for an example here.

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
                sql - one of 'STRING', 'BOOLEAN', 'INTEGER', 'REAL', or other
                    user defined type

        For example,
        {'table1':{'col1':{'primary_key':False, 'index':True, 'sql':'REAL'}}}

        INPUT:
            check -- if True, checks to make sure the database's actual structure
            matches the skeleton on record.
            EMILY - currently this throws an obnoxious exception, for debugging
            purposes. What do you think this should do, ultimately?

        EXAMPLES:
            EMILY - this would be good to do for, say, the graph database, since
            it has nontrivial structure... but to do that, it would be best if
            said database were implemented directly as a GenericSQLDatabase, so
            the example doesn't have people digging around in ext/db...

        """
        if check:
            d = construct_skeleton(self.__connection__)
            if d == self.__skeleton__:
                return d
            else:
                raise RuntimeError("BAD BAD BAD BAD BAD BAD : skeleton structure is out of whack!")
        else:
            return self.__skeleton__

    def show(self, table_name, max_field_size=20):
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

        if EMBEDDED_MODE:
            # Notebook Version
            print '<html><table bgcolor=lightgrey cellpadding=0><tr>'
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
                sql - one of 'STRING', 'BOOLEAN', 'INTEGER', 'REAL', or other
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

    TODO - mention where this database is, if it is implemented

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
        ... 'graph6':{'sql':'STRING', 'index':True, 'primary_key':True},
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
        sage: Q = SQLQuery(D, {'table_name':'simon', 'display_cols':'graph6', 'expression':['vertices','=',4]})
        sage: Q2 = SQLQuery(D, {'table_name':'simon', 'display_cols':'graph6', 'expression':['edges','=',3]})
        sage: Q = Q.intersect(Q2)
        sage: Q.run_query()
        [(u'CF', u'CF'), (u'CJ', u'CJ'), (u'CL', u'CL')]

    EMILY - explain the unicode strings, also explain why we get two copies of
    the data from the query.

    Of course, we can save the database to file:
        sage: D.save('simon.db')

    Now the database's hard link is to this file, and not the temporary db
    file. For example, let's say we open the same file with another class
    instance. We can load the file as an immutable database:
        sage: immut = GenericSQLDatabase('simon.db)
        EMILY - here is another bug, at the bottom.

    """

    def __init__(self, filename=None, skeleton=None):
        if filename is None:
            filename = tmp_filename() + '.db'
        elif (filename[-3:] != '.db'):
            raise ValueError('Please enter a valid database path (file name %s does not end in .db).'%filename)
        self.__dblocation__ = filename
        self.__connection__ = sqlite.connect(self.__dblocation__)

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
        EMILY - documentation
        """
        return self.__connection__.cursor()

    def get_connection(self):
        """
        EMILY - documentation
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
                sql - one of 'STRING', 'BOOLEAN', 'INTEGER', 'REAL', or other
                    user defined type

        EXAMPLE:
            sage: D = SQLDatabase()
            sage: table_skeleton = {
            ... 'graph6':{'sql':'STRING', 'index':True, 'primary_key':True},
            ... 'vertices':{'sql':'INTEGER'},
            ... 'edges':{'sql':'INTEGER'}
            ... }
            sage: D.create_table('simon', table_skeleton)
            sage: D.show('simon')
            edges                graph6               vertices
            ------------------------------------------------------------

        EMILY - MAKE NOTE IN DOCS THAT TO GET AN AUTO-INCREMENTING PRIMARY
        KEY, YOU MUST USE THE FULL WORD *INTEGER* INSTEAD OF *INT*, since this
        makes no sense to me.

        """
        if self.__skeleton__.has_key(table_name):
            raise ValueError("Database already has a table named %s."%table_name)

        create_statement = 'create table ' + table_name + '('
        for col in table_skeleton:
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
                sql - one of 'STRING', 'BOOLEAN', 'INTEGER', 'REAL', or other
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
        EMILY - documentation
        """
        if not self.__skeleton__.has_key(table_name):
            raise ValueError("Database has no table %s."%table_name)
        if not self.__skeleton__[table_name].has_key(col_name):
            raise ValueError("Table %s has no column %s."%(table_name,col_name))

        # TODO : bad input check?
            # check for SQL errors? something about NULL and unique==True
        # How about just a Try/Catch for more logical error messages...?
        # I'll leave this for you to decide when you are working on docstrings...
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
        EMILY - documentation
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
        EMILY - documentation

        WORD ON THE STREET IS THAT SQLITE IS RETARDED ABOUT
        *ALTER TABLE* COMMANDS... SO MEANWHILE WE ACCOMPLISH THIS
        BY CREATING A TEMPORARY TABLE.  SUGGESTIONS FOR SPEEDUP ARE
        WELCOME.  (OR JUST SEND A PATCH...)

        MAKE NOTE IN DOCS THAT TO GET AN AUTO-INCREMENTING PRIMARY
        KEY, YOU MUST USE THE FULL WORD *INTEGER* INSTEAD OF *INT*.
        """
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
        EMILY - documentation

        WORD ON THE STREET IS THAT SQLITE IS RETARDED ABOUT
        *ALTER TABLE* COMMANDS... SO MEANWHILE WE ACCOMPLISH THIS
        BY CREATING A TEMPORARY TABLE.  SUGGESTIONS FOR SPEEDUP ARE
        WELCOME.  (OR JUST SEND A PATCH...)

        Note: deosnt's delte colmsn jsut make not primiarey

        """
        if not self.__skeleton__.has_key(table_name):
            raise ValueError("Database has no table %s."%table_name)
        if not self.__skeleton__[table_name].has_key(col_name):
            raise ValueError("Table %s has no column %s."%(table_name,col_name))
        if not self.__skeleton__[table_name][index_name]['primary_key']:
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
            'create temporary table spam(%s);
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

    def add_row(self, table_name, values):
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
        if not self.__skeleton__.has_key(table_name):
            raise ValueError("Database has no table %s."%table_name)
        # values is a tuple of the right length (same as no of cols)
        if len(values) != len(self.__skeleton__[table_name]):
            raise ValueError("New row must have the same number (%d) of columns as table."%len(self.__skeleton__[table_name]))
        tup = []
        quest = "("
        for i in range (len(values)):
            tup.append(str(values[i]))
            quest += '?, '
        quest = quest.rstrip(', ') + ')'
        insert_string = 'INSERT INTO ' + table_name + ' VALUES ' + quest
        self.__connection__.execute(insert_string, tuple(tup))

##### TODO: THINK MORE ABOUT THESE WHEN LESS ILL, IE MORE BLOOD #####
#### DONE, but read through cus I left you some notes #####

    def delete_rows(self, query):
        """
        EMILY - documentation

        DOCSTRING COMMENTS:
        - Query must have no join statements (you can only delete from one table at a time)
            (This might be an interesting TODO later, ie: another function that
            creates multiple delete statements from a query with joins...  I don't think
            it's really necessary at this time though.
        - Use a query instance to modify your database.
        - Note that you cannot use a GenericSQLQuery.  (NO RESPONSIBILITY)
        - If you would like to remove all data that satisfies a query,
            enter that query as a parameter in the delete_rows function.
        - Recommended use: have some kind of primary key column that you
            use as a parameter in the query
        - Be careful, test your query first.
        - (Examples of this)
        - Also note that this database must be the database that the
            query is associated with (test that?)

        EXAMPLES:
            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_data('simon',[(0,0),(1,1),(1,2)])
            sage: p = SQLQuery(DB, {'table_name':'simon', 'display_cols':'b2', 'expression':['b2','=', 0]})
            sage: DB.delete_rows(p)
            sage: cur = DB.get_cursor()
            sage: (cur.execute('select * from simon')).fetchall()
            [(1, 1), (1, 2)]

        syntax:
        delete from table_name where blah=val
        """
        # Check query is associated with this database
        # TODO : Robert, does this work?  And do you suggest any other checking?
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

    def add_data(self, table_name, rows, entry_order=None):
        """
        INPUT:
            rows is a list of tuples that represent one row of data to add
            (types should match col types in order)
            entry_order --  an ordered list or tuple
                            overrides normal order with user defined order

        EMILY - I'm not sure we even need a try/except thing here. What I'm
        wondering is this: are there any changes to the database if SQL doesn't
        like one of the inputs? This is the only reason to do something like that.

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
        EMILY - documentation
        """
        self.__connection__.execute('vacuum')


#################################################

"""

sage: MonicPolys = SQLDatabase()
sage: MonicPolys.create_table('nums', {'n':{'sql':'INTEGER', 'index':True}})
sage: MonicPolys.add_column('nums', 'n^2', {'sql':'INTEGER', 'index':False})
Traceback (most recent call last):
...
OperationalError: unrecognized token: "^"

sage: MonicPolys.add_column('nums', 'n2', {'sql':'INTEGER', 'index':False})
Traceback (most recent call last):
...
OperationalError: unrecognized token: "^"


"""

#################################################

"""

sage: D.save('simon.db')
sage: E = SQLDatabase('simon.db')
---------------------------------------------------------------------------
<type 'exceptions.KeyError'>              Traceback (most recent call last)

/Volumes/HOME/robert/sage-2.7/<ipython console> in <module>()

/Users/robert/sage-2.7/local/lib/python2.5/site-packages/sage/databases/database.py in __init__(self, filename, skeleton)
    641
    642         # construct skeleton (from provided database)
--> 643         self.__skeleton__ = construct_skeleton(self.__connection__)
    644
    645         # add bones from new skeleton to database,

/Users/robert/sage-2.7/local/lib/python2.5/site-packages/sage/databases/database.py in construct_skeleton(connection)
    109         exe2 = cur.execute("pragma index_list(%s)"%table[0])
    110         for column in exe2.fetchall():
--> 111             skeleton[table[0]][column[1]]['index'] = True
    112     return skeleton
    113

<type 'exceptions.KeyError'>: u'sqlite_autoindex_graphs_1'

EMILY - This happens when I start a new copy of sage and try to load the same file. I'm not
sure what this means, so I don't know what to do. The D in question is the D from the tutorial
for the SQLDatabase class.

"""

################################################

"""

sage: D.save('simon.db')
sage: immut = GenericSQLDatabase('simon.db')
---------------------------------------------------------------------------
<type 'exceptions.KeyError'>              Traceback (most recent call last)

/Volumes/HOME/robert/sage-2.7/<ipython console> in <module>()

/Users/robert/sage-2.7/local/lib/python2.5/site-packages/sage/databases/database.py in __init__(self, filename)
    391         self.__connection__ = sqlite.connect(self.__dblocation__)
    392
--> 393         self.__skeleton__ = construct_skeleton(self.__connection__)
    394
    395     def copy(self):

/Users/robert/sage-2.7/local/lib/python2.5/site-packages/sage/databases/database.py in construct_skeleton(connection)
    111         exe2 = cur.execute("pragma index_list(%s)"%table[0])
    112         for column in exe2.fetchall():
--> 113             skeleton[table[0]][column[1]]['index'] = True
    114     return skeleton
    115

<type 'exceptions.KeyError'>: u'sqlite_autoindex_simon_1'


"""