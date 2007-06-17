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
        self.__query__ = ''


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

    def show(self):
        from sage.server.support import EMBEDDED_MODE
        lis = self.do_query()
        s = self.__query__[7:]
        headers = s[:s.index(' from')].split(' ')
        if len(headers) > 0:
            if headers[0] == '*':
                pass
        if EMBEDDED_MODE:  # We are in the notebook, print html table


            # split query string, ignore select and stop at join
            print '<html><table><tr>'
            print 'DATA'
            print '</html>'
        else:  # Command line, print ascii
            print lis

    def query_all(self):
        s = "select "
        pass

    def do_query(self):
        if self.__query__ == '':
            self.query_all()
        try:
            cur = self.__connection__.cursor()
            exe = cur.execute(self.__query__)
            lis = exe.fetchall()
            return lis
        except:
            raise RuntimeError('Failure to fetch query.')

    def query(self, query_string):
        self.__query__ = query_string

    def clear_query(self):
        pass
        #self.__query__ = ""