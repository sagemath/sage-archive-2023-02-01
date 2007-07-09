"""
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
            tup = str(self.__param_tuple__).rstrip(')') + ',)'
            cur = self.__database__.__connection__.cursor()
            if self.__param_tuple__ is not None:
                exe = cur.execute(self.__query_string__ + ',' + tup)
            else:
                exe = cur.execute(self.__query_string__)
            lis = exe.fetchall()
            return lis
        except:
            raise RuntimeError('Failure to fetch query.')

class SQLQuery(GenericSQLQuery):
    """
    {'table_name': 'tblname', 'display_cols': 'col1, col2, col3', 'dillhole':[col, operator, value]}
    """

    def __init__(self, database, query_dict=None):
        """
        point out strings '"value"'
        """

        if not isinstance(database, SQLDatabase):
            raise TypeError('%s is not a valid SQLDatabase'%database)

        self.__database__ = database
        self.__query_dict__ = query_dict

        #ERROR CHECKING:
        # Confirm query_dict matches skeleton:
        #      confirm tblname is in database:
        if self.__query_dict__['table_name'] not in self.__database__.__skeleton__:
                raise TypeError('%s is not a legal operator.'%operator)

        #      confirm display_cols in tblname:
        for col in self.__query_dict__[1].split(','):
            col = col.strip()
            if col not in self.__query_dict__.__skeleton__[self.__query_dict__['table_name']]:
                raise TypeError('%s column must be in %s table.'%(col,self.__query_dict__['table_name']))

        #      confirm col (from dillhole is in tblname):
        if self.__query_dict__['dillhole'][0] not in self.__query_dict__.__skeleton__[self.__query_dict__['table_name']]:
            raise TypeError('%s column must be in %s table.'%(self.__query_dict__['dillhole'][0],self.__query_dict__['table_name']))

        # confirm operator:
        verify_operator(query_dict['dillhole'][1])

        # make tuple:
        if self.__query_dict__ is not None:
            self.__param_tuple__ = (self.__query_dict__['dillhole'][2],)
        else:
            self.__param_tuple__ = None

        # make query string:
        if self.__query_dict__ is not None:
            self.__query_string__ = 'SELECT ' + self.__query_dict__['display_cols'] + \
                                    ' FROM ' + self.__query_dict__['table_name'] + \
                                    ' WHERE ' + self.__query_dict__['table_name'] + '.' + \
                                    self.__query_dict__['dillhole'][0] + ' ' + \
                                    self.__query_dict__['dillhole'][1] + ' ?'
        else:
            self.__query_string__ = None

    def intersect(self, other, join_table=None, join_dict=None):
        """
        TODO : Consider for a while and then fix the join problem
               (Apply same changes to union function)

        join_dict -- {join_table1: (corr_base_col1, col1), join_table2: (corr_base_col2, col2)}
        join_table -- base table to join on
        """

        # Check same database

        # DO SOME CHECKING - (1) if more than one table, both join args must not be None
        #                    (2) also, check in with skeleton for possible bad input
        #                    (3) and compare old query strings to confirm all previous
        #                           tables are included - (otherwise Error)
        q = SQLQuery(self.__database__)

        # inner join clause
        if join_dict is not None:
            joins = join_table
            for table in join_dict:
                joins += ' INNER JOIN ' + table + ' ON ' + join_table + \
                         '.' + table[0] + '=' + table + '.' + table[1] + ' '
            q.__query_string__ = re.sub(' FROM .* WHERE ', ' FROM' + joins + 'WHERE ', self.__query_string__)

        # concatenate display cols
        disp = self.__query_string__.split(' FROM')
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
        TODO: see intersection (look up!)
        join_dict -- {join_table1: (corr_base_col1, col1), join_table2: (corr_base_col2, col2)}
        """

        # Check same database

        # DO SOME CHECKING - (1) if more than one table, both join args must not be None
        #                    (2) also, check in with skeleton for possible bad input
        #                    (3) and compare old query strings to confirm all previous
        #                           tables are included - (otherwise Error)
        q = SQLQuery(self.__database__)

        # inner join clause
        if join_dict is not None:
            joins = join_table
            for table in join_dict:
                joins += ' INNER JOIN ' + table + ' ON ' + join_table + \
                         '.' + table[0] + '=' + table + '.' + table[1] + ' '
            q.__query_string__ = re.sub(' FROM .* WHERE ', ' FROM' + joins + 'WHERE ', self.__query_string__)

        # concatenate display cols
        disp = self.__query_string__.split(' FROM')
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
        D = Database(filename=new_loc, skeleton=copy(self.__skeleton__))
        return D

    def save(self, filename):
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

        # Silly SQLite -- we have to make a temp table to hold info...
        self.__connection__.execute('create temporary table spam(%s)'%more_attr)
        self.__connection__.execute('insert into spam select %s, %s from %s'%(original, default, table))
        self.__connection__.execute('drop table %s'%table)
        self.__connection__.execute('create table %s (%s)'%(table,more_attr))

        # Update indices in new table
        new_table_set_col_attr(self.__connection__, table, self.__skeleton__[table])

        # Now we can plop our data into the *new* table:
        self.__connection__.execute('insert into %s select %s from spam'%(table, more))
        self.__connection__.execute('drop table spam')
        self.__connection__.execute('VACUUM')

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
        self.__connection__.execute('create temporary table spam(%s)'%less_attr)
        self.__connection__.execute('insert into spam select %s from %s'%(less, table))
        self.__connection__.execute('drop table %s'%table)
        self.__connection__.execute('create table %s (%s)'%(table,less_attr))

        # Update indices in new table
        new_table_set_col_attr(self.__connection__, table, self.__skeleton__[table])

        # Now we can plop our data into the *new* table:
        self.__connection__.execute('insert into %s select %s from spam'%(table, less))
        self.__connection__.execute('drop table spam')
        self.__connection__.execute('VACUUM')

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

        # TODO : think about things like;
                # keep skeleton, don't-actually-do-anything

        self.__connection__.execute('drop table ' + table_name)

        # Update Skeleton
        self.__skeleton__.pop(table_name)

    def make_index(self, col_name, table_name, unique=False):
        if not self.__skeleton__.has_key(table_name):
            raise ValueError("Database has no table %s."%table_name)
        if not self.__skeleton__[table].has_key(col_name):
            raise ValueError("Table %s has no column %s."%(table,col_name))

        # TODO : bad input check?
            # check for SQL errors? something about NULL and unique==True
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
        self.__connection__.execute('create temporary table spam(%s)'%new_attr)
        self.__connection__.execute('insert into spam select %s from %s'%(new, table))
        self.__connection__.execute('drop table %s'%table)
        self.__connection__.execute('create table %s (%s)'%(table,new_attr))

        # Update indices in new table
        new_table_set_col_attr(self.__connection__, table, self.__skeleton__[table])

        # Now we can plop our data into the *new* table:
        self.__connection__.execute('insert into %s select %s from spam'%(table, new))
        self.__connection__.execute('drop table spam')
        self.__connection__.execute('VACUUM')

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
        self.__connection__.execute('create temporary table spam(%s)'%new_attr)
        self.__connection__.execute('insert into spam select %s from %s'%(new, table))
        self.__connection__.execute('drop table %s'%table)
        self.__connection__.execute('create table %s (%s)'%(table,new_attr))

        # Update indices in new table
        new_table_set_col_attr(self.__connection__, table, self.__skeleton__[table])

        # Now we can plop our data into the *new* table:
        self.__connection__.execute('insert into %s select %s from spam'%(table, new))
        self.__connection__.execute('drop table spam')
        self.__connection__.execute('VACUUM')

    def add_row(self, table_name, values=None):
        """
        values should be a tuple, length and order of columns in given table
        """
        if not self.__skeleton__.has_key(table_name):
            raise ValueError("Database has no table %s."%table_name)

        # EMILY!!: Read the following:
        # TODO : CHECK FOR BAD INPUT
            # -- really, sql will throw error if user is retarded enough
        # Regarding:
            # values is a tuple of the right length (same as no of cols)
        # the following will do this:
    #   if len(values) != len(self.__skeleton__[table_name]):
    #       raise ValueError("New row must have the same number (%d) of columns as table."%len(self.__skeleton__[table_name]))
        # But then why is the dafault value None?

        insert_string = 'insert into ' + table_name + ' values ' + str(values)
        cur = self.__connection__.cursor()
        exe = cur.execute(insert_string)

##### TODO: THINK MORE ABOUT THESE WHEN LESS ILL, IE MORE BLOOD #####

    def delete_rows(self, col_name, queryish_dict):
        """
        TODO :
        this is going to be really similar to the query string stuff.
        we might even want to take a query as input...  This needs more
        consideration.

        delete from table_name where blah=val
        """
        pass

    def add_data(self):
        """
        TODO :
        i.e.: from .sql file...
        still looking for this one...
        """
        pass

    def vacuum(self):
        cur = self.__connection__.cursor()
        exe = cur.execute('vacuum')


