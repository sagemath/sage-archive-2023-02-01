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
from zope.interface import implements
from twisted.internet import defer

class PublicKeyCredentialsChecker(object):
    implements(checkers.ICredentialsChecker)
    credentialInterfaces = (credentials.ISSHPrivateKey,)

    def __init__(self, pubkeydb):
        self.authorizedKeys = self.getAuthorizedKeys(pubkeydb)

    def requestAvatarId(self, credentials):
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
    credentialInterfaces = (credentials.ISSHPrivateKey,)

    def __init__(self, userdb):
        self.userdb = userdb

    def getUser(self, username):
        if not self.userdb.has_key('username'):
            return False
        return self.userdb[username]

    def requestAvatarId(self, credentials):
        user = self.getUser(credentials.username)
        if user:
            userKey = user['public_key']
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

