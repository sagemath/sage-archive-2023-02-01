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

class UserDatabase(object):
    r"""
    This class defines the UserDatabase which is used to store user
    authentication credentials and other information.

    """

    TABLENAME = 'users'
    CREATE_USER_TABLE = """CREATE TABLE %s
    (
     id integer PRIMARY KEY,
     username text NOT NULL UNIQUE,
     public_key text NOT NULL UNIQUE,
     email text,
     creation_time timestamp,
     access_time timestamp
    )
    """ % TABLENAME


    CREATION_TIME_TRIGGER = """CREATE TRIGGER creation_time
AFTER INSERT ON users
BEGIN
    UPDATE users SET creation_time = DATETIME('NOW')
    WHERE username = new.username;
END;
"""

    def __init__(self, test=False):
        self._getconf()
        self.tablename = 'users'
        if test:
            self.db_file = 'userdb_test.db'
        else:
            self.db_file = self.DB_FILE
        self.con = sqlite.connect(self.db_file,
                    detect_types=sqlite.PARSE_DECLTYPES|sqlite.PARSE_COLNAMES)
        self.con.text_factory = str

        if sql_functions.table_exists(self.con, self.tablename) is None:
            sql_functions.create_table(self.con,
                                       self.tablename,
                                       self.CREATE_USER_TABLE)
            sql_functions.add_trigger(self.con,
                                      self.CREATION_TIME_TRIGGER)
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
        except:
            print "Error reading '%s', run dsage.setup()" % conf_file
            raise
        # End reading configuration

    def get_user(self, username):
        r"""
        Returns a tuple containing the users information.

        """

        query = """SELECT * FROM users WHERE username = ?"""

        cur = self.con.cursor()
        cur.execute(query, (username,))
        result = cur.fetchone()

        return result

    def add_user(self, username, pubkey, email=None, file=True):
        r"""
        Adds a user to the database.

        """

        query = """INSERT INTO users
                   (username, public_key, email, creation_time)
                   VALUES (?, ?, ?, ?)
                """

        if file and os.path.isfile(pubkey):
            f = open(pubkey)
            type_, key = f.readlines()[0].split()[:2]
            f.close()
            if not type_ == 'ssh-rsa':
                raise TypeError
        else:
            key = pubkey

        cur = self.con.cursor()
        cur.execute(query, (username, key, email))
        self.con.commit()

    def del_user(self, username):
        r"""
        Deletes a user from the database.

        """

        query = """DELETE FROM users WHERE username = ?"""
        self.con.execute(query, username)
        self.con.commit()