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

import os
import random
import tempfile
from glob import glob

from twisted.trial import unittest
from twisted.spread import pb
from twisted.internet import reactor
from twisted.cred import portal, credentials
from twisted.conch.ssh import keys

from twisted.internet import base
base.DelayedCall.debug = True

from sage.dsage.twisted.pb import Realm
from sage.dsage.server.server import DSageServer
from sage.dsage.twisted.pb import _SSHKeyPortalRoot
from sage.dsage.twisted.pb import ClientFactory
from sage.dsage.twisted.pubkeyauth import PublicKeyCredentialsCheckerDB
from sage.dsage.database.jobdb import JobDatabaseSA as JobDatabase
from sage.dsage.database.workerdb import WorkerDatabaseSA as WorkerDatabase
from sage.dsage.database.clientdb import ClientDatabaseSA as ClientDatabase
from sage.dsage.database.db_config import init_db_sa as init_db
from sage.dsage.errors.exceptions import AuthenticationError

TEST_PUB_KEY = """ssh-rsa AAAAB3NzaC1yc2EAAAABIwAAAQEAusOUk3wZof9orc7YuKZP/wxog2uAU5BsagK4lkHgdfBc+ZR3s+Rk+k6prvuNuUXIfn2A+UkPa0xmjtQnMlqClrZXXMHhDV8iXto/vM1BopF+Ja1Y+pCK2vRRZVsZsdzL7XqyVc+kstsgKWrrguNCMIuEyc37wcsgdd1PxPmuB8Mwm3YZmNRV6yEq8Qq3IprZHfBl5S6htmwTXt4VEzvJgX1PJBLg4BauJtLxeEzYgMLY4VG3buJ2VDwlqwVPO/oVZwK3uXifXtxVx6VJO4pKUBdDSyjudPQTHxogos+8scaClx0XMh0eM7xw92j4SpA+mtzXnAKM4CqCSFH3w+/LbQ== yqiang@six
"""

TEST_PRIV_KEY = """-----BEGIN RSA PRIVATE KEY-----
MIIEoQIBAAKCAQEAusOUk3wZof9orc7YuKZP/wxog2uAU5BsagK4lkHgdfBc+ZR3
s+Rk+k6prvuNuUXIfn2A+UkPa0xmjtQnMlqClrZXXMHhDV8iXto/vM1BopF+Ja1Y
+pCK2vRRZVsZsdzL7XqyVc+kstsgKWrrguNCMIuEyc37wcsgdd1PxPmuB8Mwm3YZ
mNRV6yEq8Qq3IprZHfBl5S6htmwTXt4VEzvJgX1PJBLg4BauJtLxeEzYgMLY4VG3
buJ2VDwlqwVPO/oVZwK3uXifXtxVx6VJO4pKUBdDSyjudPQTHxogos+8scaClx0X
Mh0eM7xw92j4SpA+mtzXnAKM4CqCSFH3w+/LbQIBIwKCAQBFXpZFaJvOdM8b/F8g
A0JIyhgw0ChZjWoYvy6eNbnFZ+gE7gCTRjQidP0yXW8njvKyo6Tu4J9TvUqqFEkS
s+dcjN6ey6tcvO+CURBcEblLAtcVTwPKx/kPf1FuyhDbqcgWYMXlW8DUt8oeAyRG
juyy8f4fEf5sjUaSLaFJKYnIXc4UNTiBckcM0ffk+achcuikZKtoJ7scJVHRzS//
cz2Y/L7ma/JOASDbybwIvnZ+C+iiOrDSwCgqBWw0R48p80W3kqrPjzxMOEVtI5ON
epaM+MqWWMvpjoMBWC0ZpGl35QOLfjFjFUXzdWENfWuyjE7a6iu3S5MOBk8E21zD
OV1LAoGBAO7Df7zVUXTBNunIOlxu1fqkkAXrPUWSr9hod/kiCDM5Gtjn1519xYM7
IYtA6lipIVLVTPZlbTUd/QGxDXLloG/Ll7f0aQXz24PyF4TsVxp0cNURrsViaYEb
QCmqFqCd5fkrPGy5Qi11S+cgco2u2h6ZoyHEc4VxvMml+fv/JeG9AoGBAMg/GE56
sau4oUNKMDZGBac07uFD/IvERA/o7ul8hyuFHHRGBPIaJr7DVzDkxC0h0DscK6nl
vq2xa54yI4JH3mAPpQyCjfcT8yqcBpmDhq44udaksno7RfY7Twdk6sUMKTEbtEzk
PzRqDSIVGuYxKK4p/lPWYRkNl9AyzlDK4LNxAoGBAIGdU/jLkp53hDXEd3QBp1wt
csFicbgNzSxV9/xFrK4Xr30QJJdS55Bh7aNdwQuPA3YcBTVNAMUQR4STUHGSmOw7
UlyL/n+TAiMOZIn8pFAwlQXzqATAZSjUR2cTMNrZX5XkRV+XxNbZRnYnjqSvYHcC
8ihGEtNpoP/AgGQ6DUAHAoGAESn6xOXx+MauvJ/1gP6wBwSJgQXT0Xc5CK2Q0i8+
yTdLlPAPDW/0sUPxh9kYISAnyo1i1AxgzQ81HDAu7ejnLM4kFwPgSGDLszHx7+az
xcpZEmXjaZANT56vAKJAAkLe9ZSpDeetpWgtAuvdu/WVxckVzKv5sbC1Pbs2QXB5
qPsCgYA1GcxkJk70lyMX9W1eZ97jorqeo1+wN0IM0dAtKhB+3OCLQXEaAibtlUDL
vfFTV+I4/tOelzdrADsjxugr1RMbfBqocfaI3xJpjwNXNtqgJk8u50aFz0SlJu6y
4tizRIkoVlcZYqs9qtxEzHKryAecM8EcumKBJEJZdtdCV0O6vw==
-----END RSA PRIVATE KEY-----
"""

RANDOM_DATA =  ''.join([chr(i) for i in [random.randint(65, 123) for n in
                range(500)]])

class PublicKeyCredentialsCheckerTest(unittest.TestCase):
    username = 'tester'
    test_db = tempfile.NamedTemporaryFile()

    def setUp(self):
        Session = init_db(self.test_db.name)
        self.clientdb = ClientDatabase(Session)
        self.checker = PublicKeyCredentialsCheckerDB(self.clientdb)
        self._algorithm = 'rsa'
        pubkey = keys.Key.fromString(data=TEST_PUB_KEY)
        self._blob = keys.Key.blob(pubkey)
        self._data = RANDOM_DATA
        privkey = keys.Key.fromString(data=TEST_PRIV_KEY)
        self._signature = privkey.sign(self._data)
        self.creds = credentials.SSHPrivateKey(self.username, self._algorithm,
                                               self._blob, self._data,
                                               self._signature)
        enc_pubkey = pubkey.toString('openssh')
        self.clientdb.add_client(self.username, enc_pubkey)

    def tearDown(self):
        from sqlalchemy.orm import clear_mappers
        self.clientdb.sess.close()
        clear_mappers()
        os.remove(self.test_db.name)

    def _bad_login(self, failure):
        self.assertEquals(failure.check(AuthenticationError),
        AuthenticationError)

    def _good_login(self, result):
        self.assertEquals(result, self.username)

    def testrequestAvatarId(self):
        username = self.checker.requestAvatarId(self.creds)
        self.assertEquals(username, self.username)

        # test for bad pubkey blob
        bad_cred = credentials.SSHPrivateKey(self.username, self._algorithm,
                                             '', self._data, self._signature)
        f = self.checker.requestAvatarId(bad_cred)
        f.addErrback(self._bad_login)

        # test for bad username
        bad_cred2 = credentials.SSHPrivateKey('nonexistantuser',
                                              self._algorithm,
                                              '', self._data, '')
        f = self.checker.requestAvatarId(bad_cred2)
        f.addErrback(self._bad_login)

        # Test for a bad signature
        bad_cred3 = credentials.SSHPrivateKey(self.username, self._algorithm,
                                              self._blob, self._data, '')
        f = self.checker.requestAvatarId(bad_cred3)
        f.addErrback(self._bad_login)

        return f

    def testanonymous_login(self):
        cred = credentials.Anonymous()
        avatar_id = self.checker.requestAvatarId(cred)
        self.assertEquals(avatar_id, 'Anonymous')

    def testget_pubkey_string(self):
        from sage.dsage.twisted.pubkeyauth import get_pubkey_string
        kind, key = TEST_PUB_KEY.split()[:2]
        f = open('temp_file', 'w+b')
        f.write(TEST_PUB_KEY)
        f.close()
        self.assertEquals(key, get_pubkey_string('temp_file'))
        os.remove('temp_file')

if __name__ == 'main':
    unittest.main()
