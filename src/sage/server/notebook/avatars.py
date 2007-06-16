#####################################################################
# Copyright (C) 2007 Alex Clemesha <clemesha@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#####################################################################
from twisted.cred import portal, checkers, credentials, error as credError
from twisted.internet import protocol, defer
from zope.interface import Interface, implements
from twisted.web2 import iweb
from twisted.python import log
# import resources
#kernel connection
# import kernel


class PasswordDataBaseChecker(object):
    implements(checkers.ICredentialsChecker)
    credentialInterfaces = (credentials.IUsernamePassword,)

    def __init__(self, dbConnection):
        self.dbConnection = dbConnection

    def queryDatabase(self, result):
        if result:
            avatarId = str(result[0][0])
            return avatarId #defer.succeed(avatarId)
        else:
            return checkers.ANONYMOUS #defer.succeed(checkers.ANONYMOUS)

    def requestAvatarId(self, credentials):
        log.msg("=== requestAvatarId ===")
        username = credentials.username
        password = credentials.password
        query = "SELECT avatarId FROM users WHERE avatarId = ? AND password = ?"
        d = self.dbConnection.runQuery(query, (username, password))
        d.addCallback(self.queryDatabase)
        #d.addErrback(self._failed)
        return d

class PasswordDictChecker(object):
    implements(checkers.ICredentialsChecker)
    credentialInterfaces = (credentials.IUsernamePassword,)

    def __init__(self, passwords):
        "passwords: a dict-like object mapping usernames to passwords"
        self.passwords = passwords

    def requestAvatarId(self, credentials):
        log.msg("=== requestAvatarId ===")
        username = credentials.username
        #log.msg("un: %s, pw: %s"%(credentials.username, credentials.password))
        if self.passwords.has_key(username):
            #log.msg("password.has_key(%s)"%username)
            if credentials.password == self.passwords[username]:
                return defer.succeed(username)
            else:
                return defer.succeed(checkers.ANONYMOUS)
        else:
            return defer.succeed(checkers.ANONYMOUS)
            #return defer.fail(credError.UnauthorizedLogin("No such user"))

class LoginSystem(object):
    implements(portal.IRealm)

    def __init__(self, users):
        self.users = users #empty, stored in database right now
        # self.dbConnection = dbConnection
        self.usersResources = {} #store created resource objects
        self.kernels = {} #logged in users kernel connections.
        self.logout = lambda: None #find a good use for logout

    def requestAvatar(self, avatarId, mind, *interfaces):
        """Return a given Avatar depending on the avatarID.

        This approximatly boils down to, for a protected web site,
        that given a username (avatarId, which could just be '()' for
        an anonymous user) returned from a login page,
        (which first went through a password check in requestAvatarId)
        We serve up a given "web site" -> twisted resources, that depends
        on the avatarId, (i.e. different permissions / view depending on
        if the user is anonymous, regular, or an admin)
        """
        log.msg("=== requestAvatar ===")
        self.cookie = mind[0]
        if iweb.IResource in interfaces:
            if avatarId is checkers.ANONYMOUS: #anonymous user
                log.msg("returning AnonymousResources")
                # rsrc = resources.AnonymousRoot(self.cookie, self.dbConnection)
                from sage.server.notebook.twist import Toplevel
                rsrc = Toplevel(self.cookie)
                return (iweb.IResource, rsrc, self.logout)
            elif '@' in avatarId: #'@' in avatarId == some email address
                log.msg("returning RegularResources for %s" % avatarId)
                self._mind = mind #mind = [cookie, request.args, segments]
                self._avatarId = avatarId
                print mind[2]
                # if ('eval' or 'completer') in mind[2]:
                if ('completer' in mind[2]) or ('eval' in mind[2]):
                    self.nbid = mind[1]['nbid'][0]
                    if self.nbid in self.kernels:
                        kernelConnection = self.kernels[self.nbid]
                        print kernelConnection
                        rsrc = resources.Root(self._avatarId, self.cookie, kernelConnection, self.dbConnection)
                        return (iweb.IResource, rsrc, self.logout)
                    query = "SELECT kernel FROM notebooks WHERE notebookId = ?"
                    d = self.dbConnection.runQuery(query, (self.nbid,))
                    return d.addCallback(self.getUserResource)
                rsrc = resources.Root(avatarId, self.cookie, None, self.dbConnection)
                return (iweb.IResource, rsrc, self.logout)
        else:
            raise KeyError("None of the requested interfaces is supported")

    def getUserResource(self, result):
        ktype = str(result[0][0])
        kernelConnection = self.kernels[self.nbid] = kernel.KernelManager(ktype)
        rsrc = resources.Root(self._avatarId, self.cookie, kernelConnection, self.dbConnection)
        return (iweb.IResource, rsrc, self.logout)






