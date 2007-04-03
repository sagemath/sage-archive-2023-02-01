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

import ConfigParser
import os
import random
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

DSAGE_DIR = os.path.join(os.getenv('DOT_SAGE'), 'dsage')
# Begin reading configuration
try:
    conf_file = os.path.join(DSAGE_DIR, 'server.conf')
    config = ConfigParser.ConfigParser()
    config.read(conf_file)

    LOG_FILE = config.get('server_log', 'log_file')
    SSL = config.getint('ssl', 'ssl')
    WORKER_PORT = config.getint('server', 'worker_port')
    CLIENT_PORT = config.getint('server', 'client_port')
    PUBKEY_DATABASE = os.path.expanduser(config.get('auth',
                                                    'pubkey_database'))

    conf_file = os.path.join(DSAGE_DIR, 'client.conf')
    config = ConfigParser.ConfigParser()
    config.read(conf_file)

    LOG_FILE = config.get('log', 'log_file')
    SSL = config.getint('ssl', 'ssl')
    USERNAME = config.get('auth', 'username')
    PRIVKEY_FILE = os.path.expanduser(config.get('auth', 'privkey_file'))
    PUBKEY_FILE = os.path.expanduser(config.get('auth', 'pubkey_file'))
except:
    raise
# End reading configuration

Data =  ''.join([chr(i) for i in [random.randint(65, 123) for n in
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
        self.p.portal.registerChecker(
        PublicKeyCredentialsCheckerDB(self.clientdb))
        self.client_factory = pb.PBServerFactory(self.p)
        self.hostname = 'localhost'
        self.r = reactor.listenTCP(0, self.client_factory)
        self.port = self.r.getHost().port

        # public key authentication information
        self.username = USERNAME
        self.pubkey_file = PUBKEY_FILE
        self.privkey_file = PRIVKEY_FILE
        self.public_key_string = keys.getPublicKeyString(filename=self.pubkey_file)
        self.private_key = keys.getPrivateKeyObject(filename=self.privkey_file)
        self.public_key = keys.getPublicKeyObject(self.public_key_string)
        self.alg_name = 'rsa'
        self.blob = keys.makePublicKeyBlob(self.public_key)
        self.data = Data
        self.signature = keys.signData(self.private_key, self.data)
        self.creds = credentials.SSHPrivateKey(self.username,
                                               self.alg_name,
                                               self.blob,
                                               self.data,
                                               self.signature)
        c = ConfigParser.ConfigParser()
        c.read(os.path.join(DSAGE_DIR, 'client.conf'))
        username = c.get('auth', 'username')
        pubkey_file = c.get('auth', 'pubkey_file')
        self.clientdb.add_user(username, pubkey_file)

    def tearDown(self):
        self.connection.disconnect()
        self.jobdb._shutdown()
        files = glob('*.db*')
        for file in files:
            os.remove(file)
        return self.r.stopListening()

    def testLogin(self):
        factory = PBClientFactory()
        self.connection = reactor.connectTCP(self.hostname, self.port, factory)

        d = factory.login(self.creds, None)
        d.addCallback(self._LoginConnected)

        return d

    def _LoginConnected(self, remoteobj):
        self.assert_(isinstance(remoteobj, pb.RemoteReference))

    def testBadLogin(self):
        factory = PBClientFactory()
        self.connection = reactor.connectTCP(self.hostname, self.port, factory)

        d = factory.login(None, None)
        d.addErrback(lambda f: self.assertEquals(TypeError, f.check(TypeError)))

        return d

    def testBadLogin2(self):
        factory = PBClientFactory()
        self.connection = reactor.connectTCP(self.hostname, self.port, factory)
        bad_creds = credentials.SSHPrivateKey('this user name should not exit',
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
        self.connection = reactor.connectTCP(self.hostname, self.port, factory)
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
