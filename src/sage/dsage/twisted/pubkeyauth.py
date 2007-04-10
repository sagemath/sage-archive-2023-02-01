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

import base64

from twisted.conch import error
from twisted.conch.ssh import keys
from twisted.cred import checkers, credentials
from twisted.cred.credentials import IAnonymous
from zope.interface import implements
from twisted.internet import defer
from twisted.python import log
from twisted.spread import pb

from sage.dsage.database.clientdb import ClientDatabase
from sage.dsage.errors.exceptions import AuthenticationError

class PublicKeyCredentialsChecker(object):
    """
    This class provides authentication checking using ssh public keys.

    """

    implements(checkers.ICredentialsChecker)
    credentialInterfaces = (credentials.ISSHPrivateKey, credentials.IAnonymous)

    def __init__(self, pubkeydb):
        self.authorizedKeys = self.getAuthorizedKeys(pubkeydb)

    def requestAvatarId(self, credentials):
        if IAnonymous.providedBy(credentials):
             return 'Anonymous'

        # read the authentication table to make sure we have a fresh copy
        self.authorizedKeys = self.getAuthorizedKeys(self.file_name)

        if self.authorizedKeys.has_key(credentials.username):
            userKey = self.authorizedKeys[credentials.username]
            if not credentials.blob == base64.decodestring(userKey):
                return defer.fail(error.ConchError("Invalid key."))
            if not credentials.signature:
                return defer.fail(error.ValidPublicKey())

            pubKey = keys.getPublicKeyObject(data=credentials.blob)
            if keys.verifySignature(pubKey, credentials.signature,
                                    credentials.sigData):
                return credentials.username
            else:
                return defer.fail(error.ConchError("Invalid signature."))
        else:
            return defer.fail(error.ConchError("Invalid username."))

    def getAuthorizedKeys(self, file_name):
        self.file_name = file_name
        authorized_keys = {}
        try:
            f = open(file_name)
            for l in f:
                line = l.split(':')
                username = line[0]
                key = line[1].strip()
                authorized_keys[username] = key
        except:
            raise

        return authorized_keys

class PublicKeyCredentialsCheckerDB(object):
    implements(checkers.ICredentialsChecker)
    credentialInterfaces = (credentials.ISSHPrivateKey, credentials.IAnonymous)

    def __init__(self, clientdb):
        if not isinstance(clientdb, ClientDatabase):
            raise TypeError
        self.clientdb = clientdb

    def requestAvatarId(self, credentials):
        if IAnonymous.providedBy(credentials):
            return 'Anonymous'
        try:
            user, key = self.get_user(credentials.username)
        except TypeError:
            log.msg("Invalid username: '%s'" % credentials.username)
            return defer.fail(AuthenticationError('Login failed.'))
        if user:
            if not credentials.blob == base64.decodestring(key):
                log.msg('Invalid key.')
                return defer.fail(AuthenticationError('Login failed.'))
            if not credentials.signature:
                log.msg('No signature.')
                return defer.fail(AuthenticationError('Login failed.'))
            pub_key = keys.getPublicKeyObject(data=credentials.blob)
            if keys.verifySignature(pub_key, credentials.signature,
                                    credentials.sigData):
                # If we get to this stage, it means the user is already
                # logged in
                self.clientdb.update_login_time(credentials.username)
                return credentials.username
            else:
                log.msg('Invalid signature.')
                return defer.fail(AuthenticationError('Login failed.'))
        else:
            log.msg('Invalid username.')
            return defer.fail(AuthenticationError('Login failed.'))

    def get_user(self, username):
        return self.clientdb.get_user_and_key(username)