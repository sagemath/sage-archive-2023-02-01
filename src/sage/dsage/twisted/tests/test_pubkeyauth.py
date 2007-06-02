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
import base64
from glob import glob

from twisted.trial import unittest
from twisted.spread import pb
from twisted.internet import reactor
from twisted.cred import portal, credentials
from twisted.conch.ssh import keys

from sage.dsage.twisted.pb import Realm
from sage.dsage.server.server import DSageServer
from sage.dsage.twisted.pb import _SSHKeyPortalRoot
from sage.dsage.twisted.pb import PBClientFactory
from sage.dsage.twisted.pubkeyauth import PublicKeyCredentialsCheckerDB
from sage.dsage.database.jobdb import JobDatabaseSQLite
from sage.dsage.database.monitordb import MonitorDatabase
from sage.dsage.database.clientdb import ClientDatabase
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
    def setUp(self):
        self.jobdb = JobDatabaseSQLite(test=True)
        self.monitordb = MonitorDatabase(test=True)
        self.clientdb = ClientDatabase(test=True)
        self.dsage_server = DSageServer(self.jobdb,
                                        self.monitordb,
                                        self.clientdb,
                                        log_level=5)
        self.realm = Realm(self.dsage_server)
        self.p = _SSHKeyPortalRoot(portal.Portal(Realm(self.dsage_server)))
        self.clientdb = ClientDatabase(test=True)
        self.p.portal.registerChecker(PublicKeyCredentialsCheckerDB(
                                      self.clientdb))
        self.client_factory = pb.PBServerFactory(self.p)
        self.hostname = 'localhost'
        self.r = reactor.listenTCP(0, self.client_factory)
        self.port = self.r.getHost().port

        # public key authentication information
        self.username = 'unit_test'
        self.pubkey_file = TEST_PUB_KEY
        self.privkey_file = TEST_PRIV_KEY
        self.data = RANDOM_DATA
        self.public_key_str = keys.getPublicKeyString(data=TEST_PUB_KEY)
        self.private_key = keys.getPrivateKeyObject(
                            data=TEST_PRIV_KEY)
        self.public_key = keys.getPublicKeyObject(self.public_key_str)
        self.alg_name = 'rsa'
        self.blob = keys.makePublicKeyBlob(self.public_key)
        self.signature = keys.signData(self.private_key, self.data)
        self.creds = credentials.SSHPrivateKey(self.username,
                                               self.alg_name,
                                               self.blob,
                                               self.data,
                                               self.signature)
        pubkey = base64.encodestring(self.public_key_str).strip()
        try:
            self.clientdb.del_user(self.username)
            self.clientdb.add_user(self.username, pubkey)
        except:
            self.clientdb.add_user(self.username, pubkey)

    def tearDown(self):
        self.connection.disconnect()
        self.jobdb._shutdown()
        files = glob('*.db*')
        for file in files:
            os.remove(file)
        return self.r.stopListening()

    def testLogin(self):
        factory = PBClientFactory()
        self.connection = reactor.connectTCP(self.hostname, self.port,
                                             factory)
        d = factory.login(self.creds, None)
        d.addCallback(self._LoginConnected)

        return d

    def _LoginConnected(self, remoteobj):
        self.assert_(isinstance(remoteobj, pb.RemoteReference))

    def testBadLogin(self):
        factory = PBClientFactory()
        self.connection = reactor.connectTCP(self.hostname,
                                             self.port,
                                             factory)

        d = factory.login(None, None)
        d.addErrback(lambda f: self.assertEquals(TypeError,
                                                 f.check(TypeError)))

        return d

    def testBadLogin2(self):
        factory = PBClientFactory()
        self.connection = reactor.connectTCP(self.hostname,
                                             self.port,
                                            factory)
        bad_creds = credentials.SSHPrivateKey('bad username',
                                               self.alg_name,
                                               self.blob,
                                               self.data,
                                               self.signature)
        d = factory.login(bad_creds, None)
        d.addErrback(self._BadLoginFailure)

        return d

    def _BadLoginFailure(self, failure):
        self.assertEquals(failure.type, str(AuthenticationError))

    def testBadLogin3(self):
        factory = PBClientFactory()
        self.connection = reactor.connectTCP(self.hostname,
                                             self.port,
                                             factory)
        bad_creds = credentials.SSHPrivateKey(self.username,
                                              self.alg_name,
                                              None,
                                              self.data,
                                              self.signature)

        d = factory.login(bad_creds, None)
        d.addErrback(self._BadLoginFailure)

        return d

if __name__ == 'main':
    unittest.main()
