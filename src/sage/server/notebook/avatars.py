#####################################################################
# Copyright (C) 2007 Alex Clemesha <clemesha@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#####################################################################

SALT = 'aa'  # (it would be 4096 times more secure to not do this)

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
        return 'failed_login'
    if twist.notebook.user_is_admin(avatarId):
        return 'admin'
    return 'user'

class FailedLogin:
    def __init__(self, username):
        self.username = username

class PasswordDictChecker(object):
    implements(checkers.ICredentialsChecker)
    credentialInterfaces = (credentials.IUsernamePassword,)

    def __init__(self, passwords):
        "passwords: a dict-like object mapping usernames to passwords"
        self.passwords = passwords

    def requestAvatarId(self, credentials):
        username = credentials.username
        if self.passwords.has_key(username):
            password = self.passwords[username]
            if crypt.crypt(credentials.password, SALT) == password:
                return defer.succeed(username)
            else:
                return defer.succeed(checkers.ANONYMOUS)
        else:
            return defer.succeed(FailedLogin(username))

class PasswordFileChecker(PasswordDictChecker):
    implements(checkers.ICredentialsChecker)
    credentialInterfaces = (credentials.IUsernamePassword,)

    def __init__(self, password_file):
        """
        INPUT:
        password_file - file that contains passwords
        """

        self.password_file = password_file
        self.load_passwords()

    def load_passwords(self):
        passwords = {}
        if not os.path.exists(self.password_file):
            open(self.password_file,'w').close()
            self.add_first_admin()
        f = open(self.password_file).readlines()
        if len(f) == 0:
            self.add_first_admin()
        for line in f:
            username, password, email, account_type = line.split(':')
            password = password.strip()
            passwords[username] = password

        self.passwords = passwords

    def add_user(self, username, password, email, account_type='user'):
        self.check_username(username)
        f = open(self.password_file, 'a')
        s = '%s:%s:%s:%s\n' % (username, crypt.crypt(password, SALT), email, account_type)
        f.writelines(s)
        f.close()

    def add_first_admin(self):
        pw = "%x"%randint(2**24,2**25)
        self.add_user("admin", pw, "", "admin")

        ## TODO -- for developer convenience for now only!!!!
        log.msg("Creating accounts a and b with passwords a and b for developer convenience.")
        self.add_user("a", "a", "", "admin")
        self.add_user("b", "b", "", "user")

        log.msg("""
*************************************
     INITIALIZING USER DATABASE
*************************************
Please visit the notebook immediately
to configure the server.  Log in with
user: admin
pass: %s
*************************************"""%pw)

    def check_username(self, username):
        usernames = []
        f = open(self.password_file).readlines()
        for line in f:
            v = line.split(':')
            usernames.append(v[0])
        if username in usernames:
            raise ValueError('Username %s already exists' % username)
        else:
            return True



class LoginSystem(object):
    implements(portal.IRealm)

    def __init__(self, users):
        self.users = users #empty, stored in database right now
        # self.dbConnection = dbConnection
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
            #log.msg(avatarId)
            if avatarId is checkers.ANONYMOUS: #anonymous user
                #log.msg("returning AnonymousResources")
                rsrc = twist.AnonymousToplevel(self.cookie, avatarId)
                return (iweb.IResource, rsrc, self.logout)

            elif user_type(avatarId) == 'failed_login':
                rsrc = twist.FailedToplevel(avatarId)
                return (iweb.IResource, rsrc, self.logout)

            elif user_type(avatarId) == 'user':
                #log.msg("returning User resources for %s" % avatarId)
                self._mind = mind #mind = [cookie, request.args, segments]
                self._avatarId = avatarId
                rsrc = twist.UserToplevel(self.cookie, avatarId)
                return (iweb.IResource, rsrc, self.logout)

            elif user_type(avatarId) == 'admin':
                #log.msg("returning Admin resources for %s" % avatarId)
                self._mind = mind #mind = [cookie, request.args, segments]
                self._avatarId = avatarId
                rsrc = twist.AdminToplevel(self.cookie, avatarId)
                return (iweb.IResource, rsrc, self.logout)

        else:
            raise KeyError("None of the requested interfaces is supported")

    def getUserResource(self, result):
        ktype = str(result[0][0])
        kernelConnection = self.kernels[self.nbid] = kernel.KernelManager(ktype)
        rsrc = resources.Root(self._avatarId, self.cookie, kernelConnection, self.dbConnection)
        return (iweb.IResource, rsrc, self.logout)

