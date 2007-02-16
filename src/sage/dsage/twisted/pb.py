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

from twisted.spread import pb
from zope.interface import implements
from twisted.cred import portal, credentials
from twisted.spread.interfaces import IJellyable
from twisted.spread.pb import IPerspective, AsReferenceable
from twisted.python import log
from twisted.internet import defer

from sage.dsage.misc.hostinfo import HostInfo
from sage.dsage.server.server import DSageServer
import sage.dsage.server.client_tracker as client_tracker
import sage.dsage.server.worker_tracker as worker_tracker
from sage.dsage.errors.exceptions import BadTypeError

pb.setUnjellyableForClass(HostInfo, HostInfo)

class WorkerPBServerFactory(pb.PBServerFactory):
    def __init__(self, root):
        pb.PBServerFactory.__init__(self, root)

    def clientConnectionMade(self, broker):
        """Keeps a 3-tuple of connected workers.
        tuple[0] - the broker
        tuple[1] - the worker ip
        tuple[2] - the worker port

        """

        worker_tracker.add((broker,
                            broker.transport.getPeer().host,
                            broker.transport.getPeer().port))

    def clientConnectionLost(self, connector, reason, reconnecting=0):
        print 'Connection lost!'



class ClientPBClientFactory(pb.PBClientFactory):
    def login(self, creds, client=None):
        if not isinstance(creds, credentials.SSHPrivateKey):
            return defer.fail(TypeError())

        d = self.getRootObject()
        d.addCallback(self._cbSendUsername,
                      creds.username,
                      creds.algName,
                      creds.blob,
                      creds.sigData,
                      creds.signature,
                      client)

        return d

    def _cbSendUsername(self, root, username, alg_name, blob, sig_data,
                        signature, client):

        d = root.callRemote("login", username, alg_name, blob, sig_data,
                                signature, client)
        return d

class _SSHKeyPortalRoot(pb._PortalRoot):
    def rootObject(self, broker):
        return _SSHKeyPortalWrapper(self.portal, broker)

class _SSHKeyPortalWrapper(pb._PortalWrapper):
    def remote_login(self, username, alg_name, blob, data, signature, mind):
        pubkey_cred = credentials.SSHPrivateKey(username,
                                                alg_name,
                                                blob,
                                                data,
                                                signature)

        d = self.portal.login(pubkey_cred, mind, IPerspective)
        d.addCallback(self._loggedIn)

        return d

    def _loggedIn(self, (interface, perspective, logout)):
        if not IJellyable.providedBy(perspective):
            perspective = AsReferenceable(perspective, "perspective")
        self.broker.notifyOnDisconnect(logout)

        return perspective

class DefaultPerspective(pb.Avatar):
    r"""
    Custom implementation of pb.Avatar so we can keep track of the broker.

    """

    current_connections = 0

    def perspectiveMessageReceived(self, broker, message, args, kw):
        self.broker = broker

        return pb.Avatar.perspectiveMessageReceived(self, broker,
                                                    message, args, kw)

    def attached(self, avatar, mind):
        self.current_connections += 1
        client_tracker.add((avatar, mind))

    def detached(self, avatar, mind):
        self.current_connections -= 1
        client_tracker.remove((avatar, mind))


class UserPerspective(DefaultPerspective):
    r"""
    Defines the perspective of a regular user to the server.

    """

    def __init__(self, DSageServer, avatarID):
        self.DSageServer = DSageServer
        self.avatarID = avatarID

        log.msg('%s connected' % self.avatarID)

    def perspective_getJob(self):
        return self.DSageServer.getJob()

    def perspective_getNextJobID(self):
        job_id = self.DSageServer.getNextJobID()

        return job_id

    def perspective_getJobByID(self, jobID):
        if not isinstance(jobID, str):
            raise BadTypeError()

        job = self.DSageServer.getJobByID(jobID)

        return job

    def perspective_getJobsByAuthor(self, author, job_name, is_active):
        if not (isinstance(author, str) or
                isinstance(is_active, bool) or
                isinstance(job_name, str)):
            raise BadTypeError()

        jobs = self.DSageServer.getJobsByAuthor(author, job_name, is_active)

        return jobs

    def perspective_getJobResultByID(self, jobID):
        if not isinstance(jobID, str):
            raise BadTypeError()

        return self.DSageServer.getJobResultByID(jobID)

    def perspective_getJobOutputByID(self, jobID):
        if not isinstance(jobID, str):
            raise BadTypeError()

        return self.DSageServer.getJobOutputByID(jobID)

    def perspective_syncJob(self, jobID):
        if not isinstance(jobID, str):
            return None

        return self.DSageServer.syncJob(jobID)

    def perspective_submitJob(self, pickled_job):
        return self.DSageServer.submitJob(pickled_job)

    def perspective_jobDone(self, jobID, output, result,
                            completed, worker_info):
        if not (isinstance(jobID, str) or isinstance(completed, bool)):
            raise BadTypeError()

        return self.DSageServer.jobDone(jobID, output, result,
                                  completed, worker_info)

    def perspective_jobFailed(self, jobID):
        if not isinstance(jobID, str):
            raise BadTypeError()

        return self.DSageServer.jobFailed(jobID)

    def perspective_killJob(self, jobID, reason=None):
        if not isinstance(jobID, str):
            raise BadTypeError()

        return self.DSageServer.killJob(jobID, reason)

    def perspective_getClusterSpeed(self):
        return self.DSageServer.getClusterSpeed()

    def perspective_getWorkerList(self):
        return [x[1:] for x in self.DSageServer.getWorkerList()]

    def perspective_getClientList(self):
        return [avatar[0].avatarID for avatar in
                self.DSageServer.getClientList()]

class AdminPerspective(UserPerspective):
    r"""
    Defines the perspective of the admin.

    """

    def __init__(self, DSageServer, avatarID):
        UserPerspective.__init__(self, DSageServer, avatarID)

class Realm(object):
    implements(portal.IRealm)

    def __init__(self, DSageServer):
        self.DSageServer = DSageServer

    def requestAvatar(self, avatarID, mind, *interfaces):
        if not pb.IPerspective in interfaces:
            raise NotImplementedError, "No supported avatar interface."
        else:
            if avatarID == 'admin':
                avatar = AdminPerspective(self.DSageServer)
            else:
                avatar = UserPerspective(self.DSageServer, avatarID)
        avatar.attached(avatar, mind)
        return pb.IPerspective, avatar, lambda a=avatar:a.detached(avatar,
                                                                   mind)
