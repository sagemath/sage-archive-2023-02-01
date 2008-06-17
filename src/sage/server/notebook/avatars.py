"""nodoctest
"""
#####################################################################
# Copyright (C) 2007 Alex Clemesha <clemesha@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#####################################################################

import crypt
import os
from   random import randint

import twist
from twisted.cred import portal, checkers, credentials, error as credError
from twisted.internet import protocol, defer
from zope.interface import Interface, implements
from twisted.web2 import iweb
from twisted.python import log

def user_type(avatarId):
    if isinstance(avatarId, FailedLogin):
        if avatarId.failure_type == 'user':
            return 'invalid_user'
        elif avatarId.failure_type == 'password':
            return 'invalid_password', avatarId.username
        else:
            raise ValueError, 'invalid failure type'
    if twist.notebook.user_is_admin(avatarId):
        return 'admin'
    return 'user'

class FailedLogin:
    def __init__(self, username, failure_type = "user"):
        self.username = username
        self.failure_type = failure_type

class PasswordChecker(object):
    implements(checkers.ICredentialsChecker)
    credentialInterfaces = (credentials.IUsernamePassword,)

    def add_user(self, username, password, email, account_type='user'):
        self.check_username(username)
        U = twist.notebook.add_user(username, password, email, account_type)

    def check_username(self, username):
        usernames = twist.notebook.usernames()
        if username in usernames:
            raise ValueError('Username %s already exists' % username)
        elif username == 'pub':
            raise ValueError('"pub" is not an allowed user name')
        else:
            return True

    def requestAvatarId(self, credentials):
        username = credentials.username
        password = credentials.password
        try:
            U = twist.notebook.user(username)
        except KeyError:
            return defer.succeed(FailedLogin(username, failure_type = 'user'))

        if U.password_is(password):
            return defer.succeed(username)
        else:
            return defer.succeed(FailedLogin(username,failure_type='password'))


class ITokenCredentials(credentials.ICredentials):
    def checkToken(token):
        pass


class TokenCred(object):
    implements(ITokenCredentials)

    def __init__(self, token=None):
        self.token = token

    def checkToken(self, token):
        return token == self.token


class OneTimeTokenChecker(object):
    implements(checkers.ICredentialsChecker)
    credentialInterfaces = (ITokenCredentials,)

    is_startup = True

    def requestAvatarId(self, credentials):
        from twisted.python import log
        if self.token and credentials.checkToken(self.token):
            if self.is_startup:
                self.is_startup = False
                self.token = credentials.token = randint(0, 2**128)
            return defer.succeed('admin')
        else:
            return defer.fail("Bad token")


class LoginSystem(object):
    implements(portal.IRealm)

    def __init__(self):
        self.usersResources = {} #store created resource objects
        self.kernels = {} #logged in users kernel connections.
        self.logout = lambda: None #find a good use for logout

    def requestAvatar(self, avatarId, mind, *interfaces):
        """
        Return a given Avatar depending on the avatarID.

        This approximatly boils down to, for a protected web site,
        that given a username (avatarId, which could just be '()' for
        an anonymous user) returned from a login page,
        (which first went through a password check in requestAvatarId)
        We serve up a given "web site" -> twisted resources, that depends
        on the avatarId, (i.e. different permissions / view depending on
        if the user is anonymous, regular, or an admin)

        """
        self.cookie = mind[0]
        if iweb.IResource in interfaces:
            self._mind = mind
            self._avatarId = avatarId
            if twist.OPEN_MODE:
                rsrc = twist.AdminToplevel(self.cookie, 'admin')
                return (iweb.IResource, rsrc, self.logout)

            if avatarId is checkers.ANONYMOUS: #anonymous user
                rsrc = twist.AnonymousToplevel(self.cookie, avatarId)
                return (iweb.IResource, rsrc, self.logout)

            elif user_type(avatarId) == 'invalid_user':
                rsrc = twist.FailedToplevel(avatarId, problem='username')
                return (iweb.IResource, rsrc, self.logout)

            elif user_type(avatarId)[0] == 'invalid_password':
                rsrc = twist.FailedToplevel(avatarId, problem='password', username=user_type(avatarId)[1])
                return (iweb.IResource, rsrc, self.logout)

            elif user_type(avatarId) == 'user':
                rsrc = twist.UserToplevel(self.cookie, avatarId)
                return (iweb.IResource, rsrc, self.logout)

            elif user_type(avatarId) == 'admin':
                rsrc = twist.AdminToplevel(self.cookie, avatarId)
                return (iweb.IResource, rsrc, self.logout)

        else:
            raise KeyError("None of the requested interfaces is supported")

    def getUserResource(self, result):
        ktype = str(result[0][0])
        kernelConnection = self.kernels[self.nbid] = kernel.KernelManager(ktype)
        rsrc = resources.Root(self._avatarId, self.cookie, kernelConnection, self.dbConnection)
        return (iweb.IResource, rsrc, self.logout)

