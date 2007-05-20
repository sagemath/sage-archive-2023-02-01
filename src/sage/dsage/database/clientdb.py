##############################################################################
#
#  DSAGE: Distributed SAGE
#
#       Copyright (C) 2006, 2007 Yi Qiang <yqiang@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#
##############################################################################

import datetime
import os
import sqlite3

from twisted.python import log

from sage.dsage.twisted.pubkeyauth import get_pubkey_string
import sage.dsage.database.sql_functions as sql_functions
from sage.dsage.misc.config import get_conf

class ClientDatabase(object):
    """
    This class defines the ClientDatabase which is used to store user
    authentication credentials and other information.

    """

    TABLENAME = 'clients'

    CREATE_USER_TABLE = """CREATE TABLE %s
    (
     id integer PRIMARY KEY,
     username text NOT NULL UNIQUE,
     public_key text NOT NULL UNIQUE,
     creation_time timestamp,
     access_time timestamp,
     last_login timestamp,
     connected BOOL,
     enabled BOOL DEFAULT 1
    )
    """ % TABLENAME

    def __init__(self, test=False):
        """
        Parameters:
        test -- set to true if you would like to do testing.

        """

        self.tablename = self.TABLENAME
        if test:
            self.db_file = 'clientdb_test.db'
        else:
            self.conf = get_conf(type='clientdb')
            self.db_file = self.conf['db_file']
            if not os.path.exists(self.db_file):
                dir, file = os.path.split(self.db_file)
                if not os.path.isdir(dir):
                    os.mkdir(dir)
        self.con = sqlite3.connect(self.db_file,
                isolation_level=None,
                detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
        # Don't use this slow!
        # self.con.text_factory = sqlite3.OptimizedUnicode
        sql_functions.optimize_sqlite(self.con)
        self.con.text_factory = str
        if sql_functions.table_exists(self.con, self.tablename) is None:
            sql_functions.create_table(self.con,
                                       self.tablename,
                                       self.CREATE_USER_TABLE)
            self.con.commit()

    def _shutdown(self):
        self.con.commit()
        self.con.close()

    def get_user_and_key(self, username):
        """
        Returns a tuple containing the username and public key.

        """

        query = """SELECT username, public_key
                   FROM clients
                   WHERE username = ?"""

        cur = self.con.cursor()
        cur.execute(query, (username,))
        result = cur.fetchone()

        return result

    def get_user(self, username):
        """
        Returns a tuple containing all of the clients information.

        WARNING: ORDER OF RETURNED TUPLE MIGHT CHANGE

        """

        query = """SELECT * FROM clients WHERE username = ?"""
        cur = self.con.cursor()
        cur.execute(query, (username,))
        result = cur.fetchone()

        return result

    def add_user(self, username, pubkey):
        """
        Adds a user to the database.

        Parameters:
        username -- username
        pubkey -- public key string (must be pre-parsed already)

        """

        query = """INSERT INTO clients
                   (username, public_key, creation_time)
                   VALUES (?, ?, ?)
                """

        cur = self.con.cursor()
        cur.execute(query, (username, pubkey, datetime.datetime.now()))
        self.con.commit()

    def del_user(self, username):
        """
        Deletes a user from the database.

        """

        query = """DELETE FROM clients WHERE username = ?"""
        self.con.execute(query, (username,))
        self.con.commit()

    def set_enabled(self, username, enabled=True):
        """
        Enables/Disables a clients account.

        Parameters:
        username -- str
        enabled -- bool
        """

        query = """UPDATE clients
        SET enabled = ?
        WHERE username = ?
        """

        self.con.execute(query, (enabled, username))
        self.con.commit()

    def get_enabled(self, username):
        """
        Returns whether or not a user account is enabled.

        """

        query = """SELECT enabled FROM clients WHERE username = ?"""
        cur = self.con.cursor()
        cur.execute(query, (username,))
        result = cur.fetchone()

        return bool(result[0])

    def set_parameter(self, username, parameter, value):
        """
        Sets a particular parameter for the given username.

        """

        query = """UPDATE clients
        SET %s = ?
        WHERE username = ?
        """ % (parameter)

        self.con.execute(query, (value, username))
        self.con.commit()

    def get_parameter(self, username, parameter):
        """
        Returns a particular parameter for the given username.

        """

        query = """SELECT %s FROM clients WHERE username = ?""" % parameter
        cur = self.con.cursor()
        cur.execute(query, (username,))
        result = cur.fetchone()

        return result[0]

    def set_connected(self, username, connected=True):
        self.set_parameter(username, 'connected', connected)

    def update_login_time(self, username):
        """
        Updates the last_login time of the user.

        """

        query = """UPDATE clients
        SET last_login = ?
        WHERE username = ?
        """

        self.con.execute(query, (datetime.datetime.now(), username,))
        self.con.commit()

    def get_client_list(self):
        """
        Returns a list of clients connected.

        """

        query = """SELECT username from clients WHERE connected"""
        cur = self.con.cursor()
        cur.execute(query)

        return [result[0] for result in cur.fetchall()]