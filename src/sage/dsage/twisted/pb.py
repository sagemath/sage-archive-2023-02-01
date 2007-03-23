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
from twisted.cred.credentials import ISSHPrivateKey
from twisted.cred.credentials import Anonymous
from twisted.spread.interfaces import IJellyable
from twisted.spread.pb import IPerspective, AsReferenceable
from twisted.python import log

from sage.dsage.misc.hostinfo import HostInfo
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

        broker.notifyOnDisconnect(self.clientConnectionLost)
        worker_tracker.add((broker,
                            broker.transport.getPeer().host,
                            broker.transport.getPeer().port))

    def clientConnectionLost(self):
        for broker, host, port in worker_tracker.worker_list:
            if broker.transport.disconnected:
                worker_tracker.remove((broker, host, port))

class PBClientFactory(pb.PBClientFactory):
    r"""
    Custom implementation of the PBClientFactory that supports logging in
    with public key as well as anonymous credentials.

    """

    def login(self, creds, mind=None):
        if ISSHPrivateKey.providedBy(creds):
            d = self.getRootObject()
            d.addCallback(self._cbSendUsername,
                          creds.username,
                          creds.algName,
                          creds.blob,
                          creds.sigData,
                          creds.signature,
                          mind)

            return d
        else:
            d = self.getRootObject()
            d.addCallback(self._cbAnonymousLogin, mind)
            return d

    def _cbSendUsername(self, root, username, alg_name, blob, sig_data,
                        signature, mind):
        d = root.callRemote("login", username, alg_name, blob, sig_data,
                                signature, mind)
        return d

    def _cbAnonymousLogin(self, root, mind):
        d = root.callRemote("login_anonymous", mind)

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

    def remote_login_anonymous(self, mind):
        d = self.portal.login(Anonymous(), mind, IPerspective)
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

    def __init__(self, DSageServer, avatarID):
        self.DSageServer = DSageServer
        self.avatarID = avatarID

        log.msg('%s connected' % self.avatarID)

    def perspectiveMessageReceived(self, broker, message, args, kw):
        self.broker = broker

        return pb.Avatar.perspectiveMessageReceived(self, broker,
                                                    message, args, kw)

    def attached(self, avatar, mind):
        self.current_connections += 1
        if mind:
            self.mind = mind
            host_info = mind[1]
            host_info['ip'] = mind[0].broker.transport.getPeer().host
            host_info['port'] = mind[0].broker.transport.getPeer().port
            worker_tracker.add((mind[0].broker, host_info))
            if self.DSageServer.workerdb.get_worker(host_info['uuid']) is None:
                self.DSageServer.workerdb.add_worker(host_info)
            self.DSageServer.workerdb.set_connected(host_info['uuid'],
                                                    connected=True)
        else:
            client_tracker.add((avatar, mind))

    def detached(self, avatar, mind):
        self.current_connections -= 1
        if mind:
            # This will remove all disconnected clients, not just the one that
            # just disconnected.
            for broker, host_info in worker_tracker.worker_list:
                if broker.transport.disconnected:
                    worker_tracker.remove((broker, host_info))
            self.DSageServer.workerdb.set_connected(host_info['uuid'],
                                                    connected=False)
        else:
            client_tracker.remove((avatar, mind))

class AnonymousWorkerPerspective(DefaultPerspective):
    r"""
    Defines the perspective of an anonymous worker.

    """

    def __init__(self, DSageServer, avatarID):
        DefaultPerspective.__init__(self, DSageServer, avatarID)

    def perspective_get_job(self):
        r"""
        Returns jobs only marked as doable by anonymous workers.

        """
        return self.DSageServer.get_job(anonymous=True)

    def perspective_get_killed_jobs_list(self):
        return self.DSageServer.get_killed_jobs_list()

    def perspective_job_done(self, job_id, output,
                             result, completed, worker_info):
        if not (isinstance(job_id, str) or isinstance(completed, bool)):
            print 'Bad job_id passed to perspective_job_done'
            raise BadTypeError()

        return self.DSageServer.job_done(job_id, output, result,
                                         completed, worker_info)

class WorkerPerspective(DefaultPerspective):
    r"""
    Defines the perspective of an authenticated worker to the server.

    """
    def __init__(self, DSageServer, avatarID):
        DefaultPerspective.__init__(self, DSageServer, avatarID)

    def perspective_get_job(self):
        uuid = self.mind[1]['uuid']
        jdict = self.DSageServer.get_job(uuid=uuid)
        if jdict is not None:
            self.DSageServer.set_job_uuid(jdict['job_id'], uuid)
        return jdict

    def perspective_job_done(self, job_id, output, result,
                            completed, worker_info):
        if not (isinstance(job_id, str) or isinstance(completed, bool)):
            print 'Bad job_id passed to perspective_job_done'
            raise BadTypeError()

        return self.DSageServer.job_done(job_id, output, result,
                                  completed, worker_info)

    def perspective_job_failed(self, job_id):
        if not isinstance(job_id, str):
            print 'Bad job_id passed to perspective_job_failed'
            raise BadTypeError()

        return self.DSageServer.job_failed(job_id)

    def perspective_submit_host_info(self, hostinfo):
        if not isinstance(hostinfo, dict):
            raise BadTypeError()
        return self.DSageServer.submit_host_info(hostinfo)

    def perspective_get_killed_jobs_list(self):
        return self.DSageServer.get_killed_jobs_list()

class UserPerspective(DefaultPerspective):
    r"""
    Defines the perspective of a regular user to the server.

    """
    def __init__(self, DSageServer, avatarID):
        DefaultPerspective.__init__(self, DSageServer, avatarID)

    def perspective_get_job_by_id(self, job_id):
        if not isinstance(job_id, str):
            raise BadTypeError()
            print 'Bad job_id passed to get_job_by_id'
        log.msg('Returning job %s to %s' % (job_id, self.avatarID))
        job = self.DSageServer.get_job_by_id(job_id)

        return job

    def perspective_get_jobs_by_user_id(self, user_id):
        if not (isinstance(user_id, str)):
            print 'Bad job_id passed to perspective_get_jobs_by_user_id'
            raise BadTypeError()

        jobs = self.DSageServer.get_jobs_by_user_id(user_id)

        return jobs

    def perspective_get_job_result_by_id(self, job_id):
        if not isinstance(job_id, str):
            print 'Bad job_id passed to perspective_get_job_result_by_id'
            raise BadTypeError()

        return self.DSageServer.get_job_result_by_id(job_id)

    def perspective_get_job_output_by_id(self, job_id):
        if not isinstance(job_id, str):
            print 'Bad job_id passed to get_job_output_by_id'
            raise BadTypeError()

        return self.DSageServer.get_job_output_by_id(job_id)

    def perspective_sync_job(self, job_id):
        if not isinstance(job_id, str):
            return None

        return self.DSageServer.sync_job(job_id)

    def perspective_submit_job(self, jdict):
        return self.DSageServer.submit_job(jdict)

    def perspective_kill_job(self, job_id, reason=None):
        if not isinstance(job_id, str):
            print 'Bad job_id passed to perspective_kill_job'
            raise BadTypeError()

        return self.DSageServer.kill_job(job_id, reason)

    def perspective_get_cluster_speed(self):
        return self.DSageServer.get_cluster_speed()

    def perspective_get_worker_list(self):
        # return [x[1] for x in self.DSageServer.get_worker_list()]
        return self.DSageServer.get_worker_list()

    def perspective_get_client_list(self):
        return [avatar[0].avatarID for avatar in
                self.DSageServer.get_client_list()]

    def perspective_submit_host_info(self, hostinfo):
        if not isinstance(hostinfo, dict):
            raise BadTypeError()
        return self.DSageServer.submit_host_info(hostinfo)

    def perspective_get_killed_jobs_list(self):
        return self.DSageServer.get_killed_jobs_list()

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
                avatar = AdminPerspective(self.DSageServer, avatarID)
            elif avatarID == 'Anonymous' and mind:
                avatar = AnonymousWorkerPerspective(self.DSageServer, avatarID)
            elif mind:
                avatar = WorkerPerspective(self.DSageServer, avatarID)
            else:
                avatar = UserPerspective(self.DSageServer, avatarID)
        avatar.attached(avatar, mind)
        return pb.IPerspective, avatar, lambda a=avatar:a.detached(avatar,
                                                                   mind)