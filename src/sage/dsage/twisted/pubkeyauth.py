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

from twisted.conch.ssh import keys
from twisted.cred import checkers, credentials
from twisted.cred.credentials import IAnonymous
from zope.interface import implements
from twisted.internet import defer
from twisted.python import log

from sage.dsage.errors.exceptions import AuthenticationError

def get_pubkey_string(filename=None):
    try:
        f = open(filename)
        kind, key = f.readlines()[0].split()[:2]
        f.close()
        if not kind == 'ssh-rsa':
            raise TypeError('Invalid key type.')
    except IOError, msg:
        key = filename

    return key


class PublicKeyCredentialsCheckerDB(object):
    """
    Uses PKI to authenticate a user.
    """

    implements(checkers.ICredentialsChecker)
    credentialInterfaces = (credentials.ISSHPrivateKey,
                            credentials.IAnonymous)

    def __init__(self, clientdb):
        self.clientdb = clientdb

    def requestAvatarId(self, credentials):
        if IAnonymous.providedBy(credentials):
            return 'Anonymous'
        client = self.clientdb.get_client(credentials.username)
        if client == None:
            log.msg("Invalid username: '%s'" % credentials.username)
            return defer.fail(AuthenticationError('Login failed.'))
        try:
            pubkey = keys.Key.fromString(credentials.blob, type='blob')
        except:
            log.msg('Invalid blob for user %s' % (credentials.username))
            return defer.fail(AuthenticationError('Login failed.'))
        client_pubkey = keys.Key.fromString(client.public_key)
        if not pubkey.toString(type='openssh') == client_pubkey.toString(type='openssh'):
            log.msg('Invalid key for user %s' % (credentials.username))
            return defer.fail(AuthenticationError('Login failed.'))
        if not credentials.signature:
            log.msg('No signature for user %s ' % (credentials.username))
            return defer.fail(AuthenticationError('Login failed.'))
        if pubkey.verify(credentials.signature, credentials.sigData):
            # User is logged in at this stage
            self.clientdb.update_login_time(credentials.username)
            return credentials.username
        else:
            log.msg('Invalid signature for user %s' % (credentials.username))
            return defer.fail(AuthenticationError('Login failed.'))
