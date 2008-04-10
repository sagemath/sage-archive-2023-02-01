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
import unittest
import datetime
import os
import tempfile

from sage.dsage.database.clientdb import ClientDatabaseSQLite
from sage.dsage.database.db_config import init_db

class ClientDatabaseSQLiteTestCase(unittest.TestCase):
    """
    Test cases for DSAGE's ClientDatabaseSQLite

    """

    TEST_USERNAMES = ['yi', 'alex', 'dorian', 'robert', 'william', 'martin']
    TEST_EMAILS = ['foo@bar.com', 'bar@foo.com', '1@2.com', '2@4.com',
                   '5@6.com', '6@7.com' ]
    TEST_KEYS = ['1', '2', '3', '4', '5', '6']
    test_db = tempfile.NamedTemporaryFile()

    def setUp(self):
        db_conn = init_db(self.test_db.name)
        self.clientdb = ClientDatabaseSQLite(db_conn)

    def tearDown(self):
        query = """DELETE FROM clients"""
        cur = self.clientdb.con.cursor()
        cur.execute(query)
        self.clientdb.con.commit()
        self.clientdb._shutdown()
        os.remove(self.test_db.name)

    def testadd_client(self):
        for user, key in zip(self.TEST_USERNAMES, self.TEST_KEYS):
            self.clientdb.add_client(user, key)
            self.assert_(self.clientdb.get_client(user) is not None)

    def testdel_client(self):
        for user in self.TEST_USERNAMES:
            self.clientdb.del_client(user)
            self.assert_(self.clientdb.get_client(user) is None)

    def testset_enabled(self):
        username = self.TEST_USERNAMES[0]
        self.clientdb.add_client(username, self.TEST_KEYS[0])
        self.clientdb.set_enabled(username, enabled=True)
        self.assertEquals(self.clientdb.get_enabled(username), True)
        self.clientdb.set_enabled(username, enabled=False)
        self.assertEquals(self.clientdb.get_enabled(username), False)

    def testupdate_login_time(self):
        username = self.TEST_USERNAMES[0]
        self.clientdb.add_client(username, self.TEST_KEYS[0])
        pre_login_time = self.clientdb.get_parameter(username, 'last_login')
        self.assertEquals(pre_login_time, None)
        self.clientdb.update_login_time(username)
        post_login_time = self.clientdb.get_parameter(username, 'last_login')
        self.assert_(datetime.datetime.now() >= post_login_time)

if __name__ == '__main__':
	unittest.main()