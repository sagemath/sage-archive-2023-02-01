"""nodoctest
Keep out the bad guys.
"""

#####################################################################
# Copyright (C) 2007 Alex Clemesha <clemesha@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#####################################################################
#
# This code was inspired by the following sources:
# * The file 'guard.py' included in the Nevow Web Framework http://divmod.org/trac/wiki/DivmodNevow
# * The file 'guard.py' used in the Stiq website code http://www.stiq.it

#twisted modules
from twisted.python import log, components
from twisted.internet import task, defer
from twisted.cred.error import UnauthorizedLogin
from twisted.cred import credentials
from zope.interface import Interface, implements
from twisted.web2 import iweb
from twisted.web2 import server

#standard library
import random
import time
import md5

# SAGE's twist stuff
import twist


class Session(object):
    """A single users session, defined by the uid.
    """
    def __init__(self, uid, sessionManager):
        self.uid = uid #the cookie, a unique identifier (md5 hash)
        self.creationTime = self.lastAccessed = time.time()
        self.sessionManager = sessionManager
        self.authenticatedAs = None
        #self.setLifetime(60)

    def set_authCreds(self, creds):
        self.authenticatedAs = creds

    def get_authCreds(self):
        return self.authenticatedAs

    def get_uid(self):
        return self.uid

    def touch(self):
        self.lastAccessed = time.time()

    def expire(self):
        return defer.maybeDeferred(self.sessionManager.expiredSession, self)


class SessionsManager(object):
    """Manages all logged in users' sessions.
    """
    tickTime = 10000 #SESSION_CLEAN_FREQUENCY
    sessionLifetime = 1000000 #TRANSIENT_SESSION_LIFETIME
    sessionPersistentLifetime = 1000000 #PERSISTENT_SESSION_LIFETIME

    sessionFactory = Session

    def __init__(self):
        self.sessions = {}
        self.tick = task.LoopingCall(self._tick)

    def createSession(self):
        session = self.sessionFactory(self.createSessionID(), self)
        uid = session.get_uid()
        self.sessions[uid] = session
        #log.msg('Session %r created' % uid)
        #log.msg('All sessions %r' % self.sessions)
        if not self.tick.running and len(self.sessions) > 0:
            self.tick.start(self.tickTime)
        return session, uid

    def expiredSession(self, session):
        uid = session.get_uid()
        if uid in self.sessions:
            #log.msg('Session %r expired' % uid)
            self.sessions.pop(uid)
        if self.tick.running and len(self.sessions) == 0:
            self.tick.stop()

    def getSession(self, uid):
        session = self.sessions.get(uid)
        if session:
            session.touch()
        return session

    def createSessionID(self):
        """Generate a new session ID."""
        data = "%s_%s" % (str(random.random()) , str(time.time()))
        return md5.new(data).hexdigest()

    def _tick(self):
        """Remove expired sessions.
        """
        now = time.time()
        for uid, session in self.sessions.items():
            age = now - session.lastAccessed
            max = self.sessionLifetime
            #if session.persistent:
            #    max = self.sessionPersistentLifetime
            if age > max:
                #log.msg('Session %r expired' % uid)
                self.sessions[uid].expire()
        if not self.sessions and self.tick.running:
            self.tick.stop()

class MindManager(object):
    """Might want to use this"""
    def __init__(self, uid):
        self.uid = uid #uid is the session id (the cookie)

class MySessionWrapper(object):
    implements(iweb.IResource)

    cookieManager = None
    mindFactory = MindManager
    sessionManager = SessionsManager()

    # The interface to cred for when logging into the portal
    credInterface = iweb.IResource

    def __init__(self, portal):
        self.portal = portal

    def renderHTTP(self, request):
        """When, if ever, would this get called?
        """
        #log.msg("=== renderHTTP ===")
        d = defer.maybeDeferred(self._delegate, request, [])
        def _cb(resource, request):
            res = iweb.IResource(resource)
            return res.renderHTTP(request)
        d.addCallback(_cb, request)
        return d

    def locateChild(self, request, segments):
        """Serve Resources depending on users Authentication status.

        Inital logic occurs here to decide the
        authentication status of a given user.
        """
        if segments and segments[0] == "login":
            #log.msg("Login")
            #get the username and password in the postdata
            #the callback function needs no args because the parsing
            #of the POSTData just updates the request args.
            l = server.parsePOSTData(request)
            l.addCallback(lambda _: self.requestPasswordAuthentication(request, segments))
            return l

        if segments and segments[0] == "":
            if request.args.get('startup_token', [''])[0]:
                return self.requestPasswordAuthentication(request, segments)

        #see if the user already has a session going
        session = self.getSession(request)
        #log.msg("session: %s" % session)
        if session is None:
            #log.msg("unknown session")
            return self.requestAnonymousAuthentication(request, segments)
        else:
            if segments and segments[0] == "logout":
                #log.msg("Logout")
                return self.logout(session, request, segments)
            else:
                #log.msg("session found ... locateResource")
                creds = session.get_authCreds()
                return self.locateResource(request, segments, session, creds)

    def locateResource(self, request, segments, session, creds):
        """Locate the resource for an authenticated session.

        This method is used to actually get the resource that
        was requested after the request is passed through locateChild.

        It is in locateChild where the users session is checked.
        """
        #log.msg("=== locateResource 'myguard.py' ===")
        def _success(avatar, request, segments):
            iface, rsrc, logout = avatar
            return rsrc, segments
        #mind = self.mindFactory(request, creds)
        mind = [session.get_uid(), request.args, segments]
        d = self.portal.login(creds, mind, self.credInterface)
        d.addCallback(_success, request, segments)
        d.addErrback(self._loginFailure, request, segments, "Incorrect login.")
        return d

    def _loginFailure(self, avatar, request, segments, reason):
        return self, ()

    def _delegate(self, request, segments):
        """Identify session by the http cookie.

        If no session exists, create a new session defined
        by a uid, and then set that uid as the cookie.
        """
        cookie = get_our_cookie(request)
        if cookie in self.sessions:
            session = self.sessions[cookie]
            return self.checkLogin(request, session, segments)
        else:
            return self.sessionManager.createSession(request, segments), ()

    def requestPasswordAuthentication(self, request, segments):
        """
        Try to athenticate with a username and password from
        a web login form.  Depending on the given credentials,
        return a custom 'view' of protected resources.

        """
        #log.msg("=== requestPasswordAuthentication ===")
        creds = self.getCredentials(request)
        session, newCookie = self.sessionManager.createSession()
        mind = [newCookie, request.args, segments]
        d = self.portal.login(creds, mind, self.credInterface)
        d.addCallback(self._loginSuccess, session, creds, segments)
        return d

    def requestAnonymousAuthentication(self, request, segments):
        """
        Anonymous authentication, the user can only see
        non-protected resources.
        """
        #log.msg("=== requestAnonymousAuthentication ===")
        def _success(avatar, request, segments):
            iface, resource, logout = avatar
            return resource, segments
        creds = credentials.Anonymous() #anonymous user.
        session, newCookie = self.sessionManager.createSession()
        session.set_authCreds(creds)
        mind = [newCookie, None, None]
        d = self.portal.login(creds, mind, self.credInterface)
        d.addCallback(_success, request, segments)
        return d

    def getSession(self, request):
        """Get a user's session defined by the requests cookie.

        Pull the cookie out of the http header and
        pass it to the session manager, will handles
        find the users actual session object.
        """
        cookie = get_our_cookie(request)
        #log.msg("cookie from header: %s"%sid_cookie)
        session = self.sessionManager.getSession(cookie)
        return session

    def getCredentials(self, request):
        """Get username, password from post args
        or if they dont exist we have an anonymous
        session
        """
        if request.args.get('startup_token', [''])[0]:
            import avatars
            return avatars.TokenCred(request.args.get('startup_token', [''])[0])
        if request.headers.getHeader('cookie'):
            for C in request.headers.getHeader('cookie'):
                if C.name == 'cookie_test':
                    username = request.args.get('email', [''])[0]
                    password = request.args.get('password', [''])[0]
                else:
                    username = password = 'COOKIESDISABLED'
        else:
            username = password = 'COOKIESDISABLED'
        return credentials.UsernamePassword(username, password)

    def _loginSuccess(self, (iface, rsrc, logout), session, creds, segments):
        """
        Return the Root Page after log in success.

        Also saved the credentials that the user used to log in
        to later associate these credentials with the users session.
        """
        #log.msg("=== _loginSuccess ===")
        session.set_authCreds(creds)
        return rsrc, ()

    def _loginFailure(self, *x): #TODO
        log.msg("=== _loginFailure ===")

        log.msg(str(x))

    def incorrectLoginError(self, error, ctx, segments, loginFailure):
        pass

    def logout(self, session, request, segments):
        """Destroy the users session.

        This pops the uid that defines the users session
        out of the sessions dict.  A new anonymous session
        is immediatly created for this user.
        """
        #log.msg("=== logout ===")
        self.sessionManager.expiredSession(session)
        return self.requestAnonymousAuthentication(request, segments)





def get_our_cookie(request):
    cookies = request.headers.getHeader('cookie')
    if cookies is None:
        return None
    for C in cookies:
        if C.name == 'nb_session':
            return C.value
    return None  # not found

