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
import ConfigParser
import sqlite3 as sqlite

from twisted.python import log

import sage.dsage.database.sql_functions as sql_functions

class ClientDatabase(object):
    r"""
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
        r"""
        Parameters:
        test -- set to true if you would like to do testing.

        """

        self._getconf()
        self.tablename = self.TABLENAME
        if test:
            self.db_file = 'clientdb_test.db'
        else:
            self.db_file = self.DB_FILE
            if not os.path.exists(self.db_file):
                dir, file = os.path.split(self.db_file)
                if not os.path.isdir(dir):
                    os.mkdir(dir)
        self.con = sqlite.connect(self.db_file,
                    detect_types=sqlite.PARSE_DECLTYPES|sqlite.PARSE_COLNAMES)
        self.con.text_factory = str

        if sql_functions.table_exists(self.con, self.tablename) is None:
            sql_functions.create_table(self.con,
                                       self.tablename,
                                       self.CREATE_USER_TABLE)
            self.con.commit()

    def _shutdown(self):
        self.con.commit()
        self.con.close()

    def _getconf(self):
        self.DSAGE_DIR = os.path.join(os.getenv('DOT_SAGE'), 'dsage')
        # Begin reading configuration
        try:
            conf_file = os.path.join(self.DSAGE_DIR, 'server.conf')
            config = ConfigParser.ConfigParser()
            config.read(conf_file)

            # TODO: This needs to be changed to use db_file
            self.DB_FILE = os.path.expanduser(config.get('auth',
                                                         'pubkey_database'))
            self.LOG_FILE = config.get('db_log', 'log_file')
            self.LOG_LEVEL = config.getint('db_log', 'log_level')
        except Exception, msg:
            print msg
            print "Error reading '%s', run dsage.setup()" % conf_file
            raise
        # End reading configuration

    def get_user_and_key(self, username):
        r"""
        Returns a tuple containing the username and public key.

        """

        query = """SELECT username, public_key FROM clients WHERE username = ?"""

        cur = self.con.cursor()
        cur.execute(query, (username,))
        result = cur.fetchone()

        return result

    def get_user(self, username):
        r"""
        Returns a tuple containing all of the clients information.

        WARNING: ORDER OF RETURNED TUPLE MIGHT CHANGE

        """

        query = """SELECT * FROM clients WHERE username = ?"""
        cur = self.con.cursor()
        cur.execute(query, (username,))
        result = cur.fetchone()

        return result

    def add_user(self, username, pubkey):
        r"""
        Adds a user to the database.

        """

        query = """INSERT INTO clients
                   (username, public_key, creation_time)
                   VALUES (?, ?, ?)
                """

        try:
            f = open(pubkey)
            type_, key = f.readlines()[0].split()[:2]
            f.close()
            if not type_ == 'ssh-rsa':
                raise TypeError
        except IOError:
            key = pubkey

        cur = self.con.cursor()
        cur.execute(query, (username, key, datetime.datetime.now()))
        self.con.commit()

    def del_user(self, username):
        r"""
        Deletes a user from the database.

        """

        query = """DELETE FROM clients WHERE username = ?"""
        self.con.execute(query, (username,))
        self.con.commit()

    def set_enabled(self, username, enabled=True):
        r"""
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
        r"""
        Returns whether or not a user account is enabled.

        """

        query = """SELECT enabled FROM clients WHERE username = ?"""
        cur = self.con.cursor()
        cur.execute(query, (username,))
        result = cur.fetchone()

        return bool(result[0])

    def set_parameter(self, username, parameter, value):
        r"""
        Sets a particular parameter for the given username.

        """

        query = """UPDATE clients
        SET %s = ?
        WHERE username = ?
        """ % (parameter)

        self.con.execute(query, (value, username))
        self.con.commit()

    def get_parameter(self, username, parameter):
        r"""
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
        r"""
        Updates the last_login time of the user.

        """

        query = """UPDATE clients
        SET last_login = ?
        WHERE username = ?
        """

        self.con.execute(query, (datetime.datetime.now(), username,))
        self.con.commit()

    def get_client_list(self):
        r"""
        Returns a list of clients connected.

        """

        query = """SELECT username from clients WHERE connected"""
        cur = self.con.cursor()
        cur.execute(query)

        return [result[0] for result in cur.fetchall()]