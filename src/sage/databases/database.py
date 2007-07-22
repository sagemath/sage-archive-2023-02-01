"""
TODO : ROBERT, I am leaving you with the following (email me
                or come over if you need help):
 - renaming stuff to use a convention (as per your email)
 - renaming 'dillhole'
 - at this point, 2 try/catch blocks (i.e.: think of a fun error message)
 - docstrings and e.g.'s - (there is an example in query intersection to get you going)
 - reviewing some of my changes (look for ROBERT, TODO, EXAMPLES or DONE below)
   (also note that I threw in some examples to help you out)
 - possibly changing the add_data function to avoid max string length problem.
   (see note in function below i.e.: bottom).
 - In query show function, put in if clause to detect if user is in command line or
   notebook and decomment the notebook code.  (also in generic database show function)
 - There might be some nonetype errors throughout that should be easily
   found/fixed when making doctests.
 - Also, some query_dict (think, skeleton) checking on intersect and
   union hasn't been done yet.

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
    ABSOLUTELY NO RESPONSIBILITY FOR THIS!!!
    """

    def __init__(self, database, query_string, param_tuple=None):
        """
        TEACH PEOPLE ABOUT USING '?' AND TUPLES
        """

        if not isinstance(database, SQLDatabase):
            raise TypeError('%s is not a valid SQLDatabase'%database)

        self.__database__ = database
        self.__param_tuple__ = param_tuple
        self.__query_string__ = query_string

    def copy(self):
        return GenericSQLQuery(self.__database__, self.__query_string__, self.__param_tuple__)

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

    def run_query(self):
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
        ROBERT - Check this out!
        TODO : (Robert) - you knew of a way to tell whether SAGE
                was running in the notebook or at command line.
                Could you decomment the notebook version and throw
                in that if clause?  --ek

        EXAMPLES:
            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_data('simon',[(0,0),(1,1),(1,2)])
            sage: r = SQLQuery(DB, {'table_name':'simon', 'display_cols':'a1', 'dillhole':['b2','<=', 6]})
            sage: p = SQLQuery(DB, {'table_name':'simon', 'display_cols':'b2', 'dillhole':['b2','<=', 6]})
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

        # Notebook Version (Decomment below) - also, we
        #     could put in arguments for html keywords so the user
        #     has some control over the display
        """
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
        """

class SQLQuery(GenericSQLQuery):
    """
    query_dict := {'table_name': 'tblname', 'display_cols': 'col1, col2, col3', 'dillhole':[col, operator, value]}
    """

    def __init__(self, database, query_dict=None):
        """
        point out strings '"value"'
        """

        if not isinstance(database, SQLDatabase):
            raise TypeError('%s is not a valid SQLDatabase'%database)

        self.__database__ = database
        self.__query_dict__ = query_dict

        # TODO: ERROR CHECKING:
        # 1. Confirm query_dict matches skeleton:
        # 2. confirm tblname is in database:
        # 3. confirm display_cols in tblname:
        # 4. confirm col (from dillhole is in tblname):

        # confirm operator:
        verify_operator(query_dict['dillhole'][1])

        # make tuple:
        if self.__query_dict__ is not None:
            self.__param_tuple__ = (self.__query_dict__['dillhole'][2],)
        else:
            self.__param_tuple__ = None

        # make query string:
        if self.__query_dict__ is not None:
            self.__query_string__ = 'SELECT ' + self.__query_dict__['table_name'] + \
                                    '.' + self.__query_dict__['display_cols'] + \
                                    ' FROM ' + self.__query_dict__['table_name'] + \
                                    ' WHERE ' + self.__query_dict__['table_name'] + '.' + \
                                    self.__query_dict__['dillhole'][0] + ' ' + \
                                    self.__query_dict__['dillhole'][1] + ' ?'
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
            sage: q = SQLQuery(DB, {'table_name':'lucy', 'display_cols':'b2', 'dillhole':['a1','=',1]})
            sage: r = SQLQuery(DB, {'table_name':'simon', 'display_cols':'a1', 'dillhole':['b2','<=', 6]})
            sage: s = q.intersect(r, 'simon', {'lucy':('a1','a1')})
            sage: s.__query_string__
            'SELECT lucy.b2,simon.a1 FROM simon INNER JOIN lucy ON simon.a1=lucy.a1 WHERE ( lucy.a1 = ? ) AND ( simon.b2 <= ? )'
            sage: s.run_query()
            [(1, 1), (4, 1)]
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
        filename -- where to keep the database
    """
    def __init__(self, filename):

        if (filename[-3:] != '.db'):
            raise ValueError('Please enter a valid database path (file name %s does not end in .db).'%filename)
        self.__dblocation__ = filename
        self.__connection__ = sqlite.connect(self.__dblocation__)

        self.__skeleton__ = construct_skeleton(self.__connection__)

    def copy(self):
        """
        Returns an instance of Database with default mutable=True.
        """
        # copy .db file
        new_loc = tmp_filename() + '.db'
        os.system('cp '+ self.__dblocation__ + ' ' + new_loc)
        D = SQLDatabase(filename=new_loc, skeleton=copy(self.__skeleton__))
        return D

    def save(self, filename):
        """
        ROBERT : Read the following e.g.'s carefully.  You might want to
                change them before distributing (i.e. os.system('rm '+save_loc))
                But this should be a good intro to a bunch of functionality.
                (I mostly put it in for you).

                It covers creating a mutable db, saving it and opening it as
                an immutable db.  Then it shows modification and the new
                awesome show method.
                (This is my most recent example - the show method is probably
                    better than getting a cursor each time for the doctests so
                    I wanted to put it in since I typically did that the other
                    way -- i.e.: delete_rows method examples).
        EXAMPLES:
            sage: from sage.databases.db import DB_HOME
            sage: save_loc = DB_HOME + 'simon.db'
            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_data('simon',[(0,0),(1,1),(1,2)])
            sage: DB.save(save_loc)
            sage: immut = GenericSQLDatabase(save_loc)
            sage: immut.show('simon')
            a1                   b2
            ----------------------------------------
            0                    0
            1                    1
            1                    2
            sage: DB.show('simon')
            a1                   b2
            ----------------------------------------
            0                    0
            1                    1
            1                    2
            sage: p = SQLQuery(DB, {'table_name':'simon', 'display_cols':'b2', 'dillhole':['b2','>', 0]})
            sage: DB.delete_rows(p)
            sage: DB.show('simon')
            a1                   b2
            ----------------------------------------
            0                    0
            sage: immut.show('simon')
            a1                   b2
            ----------------------------------------
            0                    0
            1                    1
            1                    2
            sage: import os
            sage: os.system('rm ' + save_loc)
        """
        try:
            self.__connection__.execute('commit')
        except:
            # Not sure why this throws an exception - but without it,
            #       the changes are not committed so it is necessary.
            pass
        os.system('cp ' + self.__dblocation__ + ' ' + filename)

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

    def show(self, table_name, max_field_size=20):
        """
        show and entire table from the database

        TODO ROBERT:
            - same as the show function in query: throw in if clause and decomment
            - also, do you know of any quick fix for the truncated output issue in
                notebook?

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

        # Notebook Version (Decomment below) - also, we
        #     could put in arguments for html keywords so the user
        #     has some control over the display
        """
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
        """

class SQLDatabase(GenericSQLDatabase):
    """
    Dillhole and foo and piss and shit where Tom is blah and blah.

    INPUT:
        filename -- duh
        skeleton -- a triple-indexed dictionary
                outer key - table name
                    inner key - column name
                        inner inner key - one of the following:
                primary_key - boolean, whether column has been set as primary key
                index - boolean, whether column has been set as index
                sql - one of 'STRING', 'BOOLEAN', 'INTEGER', 'REAL', or other
                    user defined type


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

    def get_skeleton(self):
        # TODO (very simple):
        # debugging version of this function...
        # when done debugging, should just return instance field variable.
        return construct_skeleton(self.__connection__)

    def get_cursor(self):
        return self.__connection__.cursor()

    def get_connection(self):
        return self.__connection__

    def create_table(self, table_name, table_skeleton):
        """
        MAKE NOTE IN DOCS THAT TO GET AN AUTO-INCREMENTING PRIMARY
        KEY, YOU MUST USE THE FULL WORD *INTEGER* INSTEAD OF *INT*.

        INPUT:
            table_name -- deurrrrrrrrrrr
            table_skeleton -- a double-indexed dictionary
                outer key - column name
                    inner key - one of the following:
                primary_key - boolean, whether column has been set as primary key
                index - boolean, whether column has been set as index
                sql - one of 'STRING', 'BOOLEAN', 'INTEGER', 'REAL', or other
                    user defined type

        table_skeleton e.g.:
        {'col1': {'sql':'INTEGER', 'index':False, 'primary_key':False}, ...}

        """
        if self.__skeleton__.has_key(table_name):
            raise ValueError("Database already has a table named %s."%table_name)

        create_statement = 'create table ' + table_name + '('
        for col in table_skeleton:
            type = table_skeleton[col]['sql']
            if verify_type(type):
                if table_skeleton[col]['primary_key']:
                    create_statement += col + ' ' + type + ' primary key, '
                else:
                    create_statement += col + ' ' + type + ', '
        create_statement = create_statement.rstrip(', ') + ') '

        self.__connection__.execute(create_statement)
        new_table_set_col_attr(self.__connection__, table_name, table_skeleton)
        self.__skeleton__[table_name] = table_skeleton

    def add_column(self, table, col_name, attr_dict, default='NULL'):
        """
        Takes a while, thanks to SQLite...

        """
        # Check input:
        if not self.__skeleton__.has_key(table):
            raise ValueError("Database has no table %s."%table)
        if self.__skeleton__[table].has_key(col_name):
            raise ValueError("Table %s already has column %s."%(table,col_name))
        attr_dict = verify_column(attr_dict)

        # Get an ordered list:
        cur_list = skel_to_col_attr_list(self.__skeleton__[table])
        # Update the skeleton:
        self.__skeleton__[table][col_name] = attr_dict

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
        more_attr += col_name + ' ' + attr_dict['sql']

        # ROBERT: Look at the new fun way to do this...
        #       executescript runs a begin transaction and commit so this
        #       should speed things up for even large amounts of data

        # Silly SQLite -- we have to make a temp table to hold info...
        self.__connection__.executescript("""
            create temporary table spam(%s);
            insert into spam select %s, %s from %s;
            drop table %s;
            create table %s (%s);
            """%(more_attr, original, default, table, table, table, more_attr))

        # Update indices in new table
        new_table_set_col_attr(self.__connection__, table, self.__skeleton__[table])

        # Now we can plop our data into the *new* table:
        self.__connection__.executescript("""
            insert into %s select %s from spam;
            drop table spam;
            """%(table, more))

        self.vacuum()

    def drop_column(self, table, col_name):
        """
        Takes a while, thanks to SQLite...

        """
        # Check input:
        if not self.__skeleton__.has_key(table):
            raise ValueError("Database has no table %s."%table)
        if not self.__skeleton__[table].has_key(col_name):
            raise ValueError("Table %s has no column %s."%(table,col_name))

        # Update the skeleton:
        self.__skeleton__[table].pop(col_name)
        # Get an ordered list (without the column we're deleting):
        cur_list = skel_to_col_attr_list(self.__skeleton__[table])

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
            """%(less_attr, less, table, table, table, less_attr))
        # Update indices in new table
        new_table_set_col_attr(self.__connection__, table, self.__skeleton__[table])

        # Now we can plop our data into the *new* table:
        self.__connection__.executescript("""
            insert into %s select %s from spam;
            drop table spam;
            """%(table, less))

        self.vacuum()

    def rename_table(self, table, new_name):
        # Check input:
        if not self.__skeleton__.has_key(table):
            raise ValueError("Database has no table %s."%table)
        if self.__skeleton__.has_key(new_name):
            raise ValueError("Database already has table %s."%new_name)

        self.__connection__.execute('alter table %s rename to %s'%(table, new_name))

        # Update skeleton:
        self.__skeleton__[new_name] = self.__skeleton__[table]
        self.__skeleton__.pop(table)

    def drop_table(self, table_name):
        if not self.__skeleton__.has_key(table_name):
            raise ValueError("Database has no table %s."%table_name)

        # (This TODO is DONE): think about things like;
                # keep skeleton, don't-actually-do-anything
        # DONE : handled by creating new function
                # drop_data_from_table (see below)

        self.__connection__.execute('drop table ' + table_name)

        # Update Skeleton
        self.__skeleton__.pop(table_name)

    def drop_data_from_table(self, table_name):
        """
        skeleton stays same but all data disappears...
        """
        if not self.__skeleton__.has_key(table_name):
            raise ValueError("Database has no table %s."%table_name)
        self.__connection__.execute('delete from ' + table_name)

    def make_index(self, col_name, table_name, unique=False):
        if not self.__skeleton__.has_key(table_name):
            raise ValueError("Database has no table %s."%table_name)
        if not self.__skeleton__[table].has_key(col_name):
            raise ValueError("Table %s has no column %s."%(table,col_name))

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

    def make_primary_key(self, table, col_name):
        """
        WORD ON THE STREET IS THAT SQLITE IS RETARDED ABOUT
        *ALTER TABLE* COMMANDS... SO MEANWHILE WE ACCOMPLISH THIS
        BY CREATING A TEMPORARY TABLE.  SUGGESTIONS FOR SPEEDUP ARE
        WELCOME.  (OR JUST SEND A PATCH...)

        MAKE NOTE IN DOCS THAT TO GET AN AUTO-INCREMENTING PRIMARY
        KEY, YOU MUST USE THE FULL WORD *INTEGER* INSTEAD OF *INT*.
        """
        if not self.__skeleton__.has_key(table):
            raise ValueError("Database has no table %s."%table)
        if not self.__skeleton__[table].has_key(col_name):
            raise ValueError("Table %s has no column %s."%(table,col_name))

        # Update the skeleton:
        self.__skeleton__[table][col_name]['primary_key'] = True
        # Get an ordered list (with the primary key info updated):
        cur_list = skel_to_col_attr_list(self.__skeleton__[table])

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
            """%(new_attr, new, table, table, table,new_attr))

        # Update indices in new table
        new_table_set_col_attr(self.__connection__, table, self.__skeleton__[table])

        # Now we can plop our data into the *new* table:
        self.__connection__.executescript("""
            insert into %s select %s from spam;
            drop table spam;
            """%(table, new))
        self.vacuum()

    def drop_primary_key(self, table, col_name):
        """
        WORD ON THE STREET IS THAT SQLITE IS RETARDED ABOUT
        *ALTER TABLE* COMMANDS... SO MEANWHILE WE ACCOMPLISH THIS
        BY CREATING A TEMPORARY TABLE.  SUGGESTIONS FOR SPEEDUP ARE
        WELCOME.  (OR JUST SEND A PATCH...)

        Note: deosnt's delte colmsn jsut make not primiarey

        """
        if not self.__skeleton__.has_key(table):
            raise ValueError("Database has no table %s."%table)
        if not self.__skeleton__[table].has_key(col_name):
            raise ValueError("Table %s has no column %s."%(table,col_name))
        if not self.__skeleton__[table_name][index_name]['primary_key']:
            return # silently

        # Update the skeleton:
        self.__skeleton__[table][col_name]['primary_key'] = False
        # Get an ordered list (with the primary key info updated):
        cur_list = skel_to_col_attr_list(self.__skeleton__[table])

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
            """%(new_attr, new, table, table, table, new_attr))

        # Update indices in new table
        new_table_set_col_attr(self.__connection__, table, self.__skeleton__[table])

        # Now we can plop our data into the *new* table:
        self.__connection__.executescript("""
            insert into %s select %s from spam;
            drop table spam;
            """%(table, new))
        self.vacuum()

    def add_row(self, table_name, values):
        """
        values should be a tuple, length and order of columns in given table

        EXAMPLES:
            sage: DB = SQLDatabase()sage: DB.create_table('simon',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_row('simon',(0,1))
            sage: cur = DB.get_cursor()
            sage: (cur.execute('select * from simon')).fetchall()
            [(0, 1)]
        """
        if not self.__skeleton__.has_key(table_name):
            raise ValueError("Database has no table %s."%table_name)

        # ROBERT!!: Read the following:  (I decommented some, and changed the default
        #           value because it should not have been None)
        # TODO : CHECK FOR BAD INPUT
        #  -- I would throw in a try/catch block here...
            # -- really, sql will throw error if user is retarded enough

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
        ROBERT : I haven't tested this... if if doesn't work, then it is probably the
            param_tuple thing...  send me an email if it's not working

        EXAMPLES:
            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool','primary_key':False}, 'b2':{'sql':'int', 'primary_key':False}})
            sage: DB.add_data('simon',[(0,0),(1,1),(1,2)])
            sage: p = SQLQuery(DB, {'table_name':'simon', 'display_cols':'b2', 'dillhole':['b2','=', 0]})
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
        if (query.__query_string__).__contains__(' JOIN '):
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

        ROBERT: should we do a try/catch
            or check each type?  Probably, try/catch (faster than checking all).
            That would make three then.  TODO.
        ROBERT: I'm a little worried about max string length in python.  Could you
            look into a way around that particular error?  i.e.: writing to and then
            reading from a temporary file or something... (instead of string *script*)

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
        self.__connection__.execute('vacuum')


