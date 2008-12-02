"""
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
    """
    Return the type of user specified by the given avatarId, which is
    either a string or an instance of the FailedLogin class.

    INPUT:
        avatarId -- string or FailedLogin instance
    OUTPUT:
        string -- 'invalid_user', 'admin', 'user'

    EXAMPLES:
        sage: import sage.server.notebook.twist
        sage: import sage.server.notebook.avatars as avatars
        sage: avatars.user_type(avatars.FailedLogin('fake'))
        'invalid_user'
        sage: avatars.user_type('_sage_')
        'invalid_user'
        sage: avatars.user_type('pub')
        'invalid_user'
        sage: avatars.user_type('guest')
        'invalid_user'
        sage: avatars.user_type('admin')
        Traceback (most recent call last):
        ...
        AttributeError: 'NoneType' object has no attribute 'user_is_admin'
        sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
        sage: nb.create_default_users('password')
        Creating default users.
        sage: sage.server.notebook.twist.notebook = nb
        sage: avatars.user_type('admin')
        'admin'
        sage: nb.add_user('bob', 'an**d', 'bob@gmail.com', force=True)
        sage: avatars.user_type('bob')
        'user'
    """
    if isinstance(avatarId, FailedLogin):
        if avatarId.failure_type == 'user':
            return 'invalid_user'
        elif avatarId.failure_type == 'password':
            return 'invalid_password', avatarId.username
        elif avatarId.failure_type == 'cookies':
            return 'cookies_disabled'
        else:
            raise ValueError, 'invalid failure type'

    # It is critically important that it be impossible to login as the
    # pub, _sage_, or guest users.  This _sage_ user is a fake user that is used
    # internally by the notebook for the doc browser and other tasks.
    if avatarId in ['_sage_', 'guest', 'pub']:
        return 'invalid_user'

    # This only works once the notebook object in the twist
    # module has been initialized, which happens only once
    # the notebook starts running.
    if isinstance(avatarId, str) and twist.notebook.user_is_admin(avatarId):
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
        if username == 'COOKIESDISABLED':
            return defer.succeed(FailedLogin(username, failure_type = 'cookies'))

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

        This approximately boils down to, for a protected web site,
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

            T = user_type(avatarId)

            if T == 'invalid_user':
                rsrc = twist.FailedToplevel(avatarId, problem='username')
                return (iweb.IResource, rsrc, self.logout)

            elif T[0] == 'invalid_password':
                rsrc = twist.FailedToplevel(avatarId, problem='password', username=user_type(avatarId)[1])
                return (iweb.IResource, rsrc, self.logout)

            elif T == 'cookies_disabled':
                rsrc = twist.FailedToplevel(avatarId, problem='cookies_disabled')
                return (iweb.IResource, rsrc, self.logout)

            elif T == 'user':
                rsrc = twist.UserToplevel(self.cookie, avatarId)
                return (iweb.IResource, rsrc, self.logout)

            elif T == 'admin':
                rsrc = twist.AdminToplevel(self.cookie, avatarId)
                return (iweb.IResource, rsrc, self.logout)

        else:
            raise KeyError("None of the requested interfaces is supported")

    def getUserResource(self, result):
        ktype = str(result[0][0])
        kernelConnection = self.kernels[self.nbid] = kernel.KernelManager(ktype)
        rsrc = resources.Root(self._avatarId, self.cookie, kernelConnection, self.dbConnection)
        return (iweb.IResource, rsrc, self.logout)

