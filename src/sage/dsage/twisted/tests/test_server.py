############################################################################
#
#   DSAGE: Distributed SAGE
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
############################################################################

#import unittest
#import os

#class ServerTestCase(unittest.TestCase):
#    def testgetAuthorizedKeys(self):
#        cur = os.getcwd()
#        os.chdir('../dsage/tests/data')
#        authorized_keys = getAuthorizedKeys('authorized_keys.db')
#        for username, key in authorized_keys.iteritems():
#            self.assert_(isinstance(username, str))
#            self.assert_(isinstance(key, str))
#        os.chdir(cur)
#
