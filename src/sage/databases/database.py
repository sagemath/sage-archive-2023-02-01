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
from sqlite3 import dbapi2 as sqlite
import os
from sage.misc.misc import tmp_filename
from sage.structure.sage_object import SageObject

# define regexp?

# define other fcns?

class ImmutableDatabase(SageObject):

    ## __dblocation__, __connection__, __query__, __skeleton__, __auto_key__

    def __init__(self, filename):
        """
        Immutable Database.

        INPUT:
            filename -- where to keep the database


        """
        if (filename[-3:] != '.db'):
            raise ValueError('Please enter a valid database path (file name %s does not end in .db).'%filename)
        self.__dblocation__ = filename
        self.__connection__ = sqlite.connect(self.__dblocation__)

        # construct skeleton
        self.__skeleton__ = {}
        cur = self.__connection__.cursor()
        exe = cur.execute("select name from sqlite_master where type='table'")
        for table in exe.fetchall():
            self.__skeleton__[table[0]] = {}
            exe1 = cur.execute("pragma table_info(%s)"%table[0])
            for column in exe1.fetchall():
                self.__skeleton__[table[0]][column[1]] = {'sql':column[2]}



#        exe = cur.execute(""


    def copy():
        """
        Returns an instance of Database with default mutable=True.
        Deep copy of skeleton, auto_key and query with new dblocation.
        """
        # copy .db file
        new_loc = tmp_filename() + '.db'
        os.system('cp '+ self.__dblocation__ + ' ' + new_loc)
        D = Database(new_loc, copy(self.__skeleton__), copy(self.__auto_key__))
        return D

    def save():
        pass

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

    def show():
        pass

    def query():
        pass

    def clear_query():
        pass