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

from sage.dsage.database.userdb import UserDatabase

class UserDatabaseTestCase(unittest.TestCase):
    r"""
    Test cases for DSAGE's UserDatabase

    """

    TEST_USERNAMES = ['yi', 'alex', 'dorian', 'robert', 'william', 'martin']
    TEST_EMAILS = ['foo@bar.com', 'bar@foo.com', '1@2.com', '2@4.com', '5@6.com', '6@7.com' ]
    TEST_KEYS = ['1', '2', '3', '4', '5', '6']

    def setUp(self):
        self.userdb = UserDatabase(test=True)

    def tearDown(self):
        self.userdb._shutdown()

    def testadd_user(self):
        for user, key, email in zip(self.TEST_USERNAMES, self.TEST_KEYS, self.TEST_EMAILS):
            self.userdb.add_user(user, key, email)

if __name__ == '__main__':
	unittest.main()