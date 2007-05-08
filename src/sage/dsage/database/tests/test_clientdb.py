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

from sage.dsage.database.clientdb import ClientDatabase

class ClientDatabaseTestCase(unittest.TestCase):
    """
    Test cases for DSAGE's ClientDatabase

    """

    TEST_USERNAMES = ['yi', 'alex', 'dorian', 'robert', 'william', 'martin']
    TEST_EMAILS = ['foo@bar.com', 'bar@foo.com', '1@2.com', '2@4.com', '5@6.com', '6@7.com' ]
    TEST_KEYS = ['1', '2', '3', '4', '5', '6']

    def setUp(self):
        self.clientdb = ClientDatabase(test=True)

    def tearDown(self):
        query = """DELETE FROM clients"""
        cur = self.clientdb.con.cursor()
        cur.execute(query)
        self.clientdb.con.commit()
        self.clientdb._shutdown()

    def testadd_user(self):
        for user, key in zip(self.TEST_USERNAMES, self.TEST_KEYS):
            self.clientdb.add_user(user, key)
            self.assert_(self.clientdb.get_user(user) is not None)

    def testdel_user(self):
        for user in self.TEST_USERNAMES:
            self.clientdb.del_user(user)
            self.assert_(self.clientdb.get_user(user) is None)

    def testset_enabled(self):
        username = self.TEST_USERNAMES[0]
        self.clientdb.add_user(username, self.TEST_KEYS[0])
        self.clientdb.set_enabled(username, enabled=True)
        self.assertEquals(self.clientdb.get_enabled(username), True)
        self.clientdb.set_enabled(username, enabled=False)
        self.assertEquals(self.clientdb.get_enabled(username), False)

    def testupdate_login_time(self):
        username = self.TEST_USERNAMES[0]
        self.clientdb.add_user(username, self.TEST_KEYS[0])
        pre_login_time = self.clientdb.get_parameter(username, 'last_login')
        self.assertEquals(pre_login_time, None)
        self.clientdb.update_login_time(username)
        post_login_time = self.clientdb.get_parameter(username, 'last_login')
        self.assert_(datetime.datetime.now() >= post_login_time)

if __name__ == '__main__':
	unittest.main()