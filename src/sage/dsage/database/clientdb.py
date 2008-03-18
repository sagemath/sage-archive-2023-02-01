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

from sage.dsage.misc.constants import SERVER_LOG
from sage.dsage.database.client import Client

class ClientDatabaseSA(object):
    def __init__(self, Session):
        self.sess = Session()

    def _add_test_client(self):
        print 'Adding testing client...'
        self.add_client('tester', '')

    def get_client(self, username):
        c = self.sess.query(Client).filter_by(username=username).first()

        return c

    def add_client(self, username, public_key):
        c = Client(username, public_key)
        self.sess.save(c)
        self.sess.commit()

    def del_client(self, username):
        c = self.sess.query(Client).filter_by(username=username).first()
        self.sess.delete(c)
        self.sess.commit()

    def get_client_list(self):
        clients = self.sess.query(Client).all()

        return clients

    def update_login_time(self, username):
        c = self.sess.query(Client).filter_by(username=username).first()
        c.last_login = datetime.datetime.now()
        self.sess.save_or_update(c)
        self.sess.commit()

    def set_connected(self, username, connected):
        c = self.sess.query(Client).filter_by(username=username).first()
        c.connected = connected
        self.sess.save_or_update(c)
        self.sess.commit()

class ClientDatabaseSQLite(object):
    """
    This class defines the ClientDatabase which is used to store user
    authentication credentials and other information.

    """

    def __init__(self, db_conn, log_file=SERVER_LOG, log_level=0):
        """
        Parameters:
        test -- set to true if you would like to do testing.

        """

        self.log_level = log_level
        self.log_file = log_file
        self.tablename = 'clients'
        self.con = db_conn

    def _shutdown(self):
        self.con.commit()

    def get_client_key(self, username):
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

    def get_client(self, username):
        """
        Returns a tuple containing all of the clients information.

        WARNING: ORDER OF RETURNED TUPLE MIGHT CHANGE

        """

        query = """SELECT * FROM clients WHERE username = ?"""
        cur = self.con.cursor()
        cur.execute(query, (username,))
        result = cur.fetchone()

        return result

    def add_client(self, username, pubkey):
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

    def del_client(self, username):
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

    def get_client_list(self, connected=True):
        """
        Returns a list of clients connected.

        """

        if connected:
            query = """SELECT * from clients WHERE connected"""
        else:
            query = """SELECT * from clients"""
        cur = self.con.cursor()
        cur.execute(query)
        result = cur.fetchall()

        return [self.create_dict(r, cur.description) for r in result]

    def create_dict(self, tuple_, row_description):
        columns = [desc[0] for desc in row_description]
        d = dict(zip(columns, tuple_))

        return d
