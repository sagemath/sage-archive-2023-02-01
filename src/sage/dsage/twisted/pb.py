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
from twisted.spread import banana
banana.SIZE_LIMIT = 100*1024*1024 # 100 MegaBytes

from zope.interface import implements
from twisted.cred import portal, credentials
from twisted.cred.credentials import ISSHPrivateKey
from twisted.cred.credentials import Anonymous
from twisted.spread.interfaces import IJellyable
from twisted.spread.pb import IPerspective, AsReferenceable
from twisted.python import log

from sage.dsage.misc.hostinfo import HostInfo
import sage.dsage.server.worker_tracker as worker_tracker
from sage.dsage.errors.exceptions import BadTypeError, BadJobError

pb.setUnjellyableForClass(HostInfo, HostInfo)

class WorkerPBServerFactory(pb.PBServerFactory):
    """
    This factory serves workers requests.

    """

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
    """
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
    """
    Custom implementation of pb.Avatar.

    """

    def __init__(self, DSageServer, avatarID):
        self.DSageServer = DSageServer
        self.avatarID = avatarID
        self.connections = 0

    def perspectiveMessageReceived(self, broker, message, args, kw):
        self.broker = broker

        return pb.Avatar.perspectiveMessageReceived(self, broker,
                                                    message, args, kw)

    def attached(self, avatar, mind):
        self.connections += 1
        if isinstance(mind, tuple):
            self.mind = mind
            self.host_info = mind[1]
            self.host_info['ip'] = mind[0].broker.transport.getPeer().host
            self.host_info['port'] = mind[0].broker.transport.getPeer().port
            uuid = self.host_info['uuid']
            if self.DSageServer.monitordb.get_monitor(uuid) is None:
                self.DSageServer.monitordb.add_monitor(self.host_info)
            self.DSageServer.monitordb.set_connected(uuid, connected=True)
        else:
            self.DSageServer.clientdb.update_login_time(self.avatarID)
            self.DSageServer.clientdb.set_connected(self.avatarID,
                                                    connected=True)

    def detached(self, avatar, mind):
        self.connections -= 1
        log.msg('%s disconnected' % (self.avatarID))
        if isinstance(mind, tuple):
            self.DSageServer.monitordb.set_connected(self.host_info['uuid'],
                                                     connected=False)
        else:
            self.DSageServer.clientdb.set_connected(self.avatarID,
                                                    connected=False)

class AnonymousMonitorPerspective(DefaultPerspective):
    """
    Defines the perspective of an anonymous worker.

    """

    def __init__(self, DSageServer, avatarID):
        DefaultPerspective.__init__(self, DSageServer, avatarID)

    def attached(self, avatar, mind):
        self.connections += 1
        self.mind = mind
        self.host_info = mind[1]
        self.host_info['ip'] = mind[0].broker.transport.getPeer().host
        self.host_info['port'] = mind[0].broker.transport.getPeer().port
        uuid = self.host_info['uuid']
        if self.DSageServer.monitordb.get_monitor(uuid) is None:
            self.DSageServer.monitordb.add_monitor(self.host_info)
        self.DSageServer.monitordb.set_connected(uuid, connected=True)
        self.DSageServer.monitordb.set_anonymous(uuid, anonymous=True)

    def detached(self, avatar, mind):
        self.connections -= 1
        log.msg('%s disconnected' % (self.avatarID))
        self.DSageServer.monitordb.set_connected(self.host_info['uuid'],
                                                 connected=False)

    def perspective_get_job(self):
        """
        Returns jobs only marked as doable by anonymous workers.

        """

        uuid = self.mind[1]['uuid']
        jdict = self.DSageServer.get_job(anonymous=True)
        if jdict is not None:
            self.DSageServer.set_job_uuid(jdict['job_id'], uuid)
            self.DSageServer.set_busy(uuid, busy=True)
        else:
            self.DSageServer.set_busy(uuid, busy=False)

        return jdict

    def perspective_get_killed_jobs_list(self):
        return self.DSageServer.get_killed_jobs_list()


    def perspective_job_failed(self, job_id, traceback):
        if not isinstance(job_id, str):
            log.msg('Bad job_id %s' % (job_id))
            raise BadTypeError()

        uuid = self.mind[1]['uuid']
        self.DSageServer.set_busy(uuid, busy=False)

        return self.DSageServer.job_failed(job_id, traceback)

    def perspective_job_done(self, job_id, output, result, completed):
        if not (isinstance(job_id, str) or isinstance(completed, bool)):
            log.msg('Bad job_id passed to perspective_job_done')
            log.msg('job_id: %s' % (job_id))
            log.msg('output: %s' % (output))
            log.msg('completed: %s' % (completed))
            raise BadTypeError()
        if completed:
            uuid = self.mind[1]['uuid']
            self.DSageServer.set_busy(uuid, busy=False)

        return self.DSageServer.job_done(job_id, output, result, completed)

    def perspective_submit_host_info(self, hostinfo):
        if not isinstance(hostinfo, dict):
            raise BadTypeError()

        return self.DSageServer.submit_host_info(hostinfo)

class MonitorPerspective(AnonymousMonitorPerspective):
    """
    Defines the perspective of an authenticated worker to the server.

    """
    def __init__(self, DSageServer, avatarID):
        DefaultPerspective.__init__(self, DSageServer, avatarID)

    def attached(self, avatar, mind):
        self.connections += 1
        self.mind = mind
        self.host_info = mind[1]
        self.host_info['ip'] = mind[0].broker.transport.getPeer().host
        self.host_info['port'] = mind[0].broker.transport.getPeer().port
        uuid = self.host_info['uuid']
        if self.DSageServer.monitordb.get_monitor(uuid) is None:
            self.DSageServer.monitordb.add_monitor(self.host_info)
        self.DSageServer.monitordb.set_connected(uuid, connected=True)
        self.DSageServer.monitordb.set_anonymous(uuid, anonymous=False)

    def perspective_get_job(self):
        """
        Returns jobs to authenticated workers.

        """

        try:
            uuid = self.mind[1]['uuid']
        except Exception, msg:
            raise ValueError("Could not match a uuid to the monitor.")

        jdict = self.DSageServer.get_job(anonymous=False)
        if jdict is not None:
            self.DSageServer.set_job_uuid(jdict['job_id'], uuid)
            self.DSageServer.set_busy(uuid, busy=True)
        else:
            self.DSageServer.set_busy(uuid, busy=False)

        return jdict

class UserPerspective(DefaultPerspective):
    """
    Defines the perspective of a regular user to the server.

    """

    def __init__(self, DSageServer, avatarID):
        DefaultPerspective.__init__(self, DSageServer, avatarID)

    def perspective_get_job_by_id(self, job_id):
        if not isinstance(job_id, str):
            log.msg('Bad job_id [%s] passed to get_job_by_id' % (job_id))
            raise BadTypeError()
        # log.msg('Returning job %s to %s' % (job_id, self.avatarID))
        job = self.DSageServer.get_job_by_id(job_id)

        return job

    def perspective_get_jobs_by_username(self, username, active=True):
        if not (isinstance(username, str)):
            log.msg('Bad username [%s] passed to perspective_get_jobs_by_username' % (username))
            raise BadTypeError()

        jobs = self.DSageServer.get_jobs_by_username(username, active)

        return jobs

    def perspective_get_job_result_by_id(self, job_id):
        if not isinstance(job_id, str):
            log.msg('Bad job_id [%s] passed to perspective_get_job_result_by_id' % (job_id))
            raise BadTypeError()

        return self.DSageServer.get_job_result_by_id(job_id)

    def perspective_get_job_output_by_id(self, job_id):
        if not isinstance(job_id, str):
            log.msg('Bad job_id [%s] passed to get_job_output_by_id' % (job_id))
            raise BadTypeError()

        return self.DSageServer.get_job_output_by_id(job_id)

    def perspective_sync_job(self, job_id):
        if not isinstance(job_id, str):
            return None

        return self.DSageServer.sync_job(job_id)

    def perspective_submit_job(self, jdict):
        if jdict is None:
            raise BadJobError()
        if jdict['username'] != self.avatarID:
            log.msg('username does not match credentials')
            log.msg('claim: %s' % jdict['username'])
            log.msg('actual: %s' % self.avatarID)
            raise BadJobError()

        return self.DSageServer.submit_job(jdict)

    def perspective_kill_job(self, job_id, reason=None):
        if not isinstance(job_id, str):
            log.msg('Bad job_id [%s] passed to perspective_kill_job' % job_id)
            raise BadTypeError()

        return self.DSageServer.kill_job(job_id, reason)

    def perspective_get_cluster_speed(self):
        return self.DSageServer.get_cluster_speed()

    def perspective_get_monitor_list(self):
        # return [x[1] for x in self.DSageServer.get_worker_list()]
        return self.DSageServer.get_monitor_list()

    def perspective_get_client_list(self):
        return self.DSageServer.get_client_list()

    def perspective_get_worker_count(self):
        return self.DSageServer.get_worker_count()

    def perspective_submit_host_info(self, hostinfo):
        if not isinstance(hostinfo, dict):
            raise BadTypeError()
        return self.DSageServer.submit_host_info(hostinfo)

    def perspective_get_killed_jobs_list(self):
        return self.DSageServer.get_killed_jobs_list()

class AdminPerspective(UserPerspective):
    """
    Defines the perspective of the admin.

    """

    def __init__(self, DSageServer, avatarID):
        UserPerspective.__init__(self, DSageServer, avatarID)

class Realm(object):
    implements(portal.IRealm)

    def __init__(self, DSageServer):
        self.DSageServer = DSageServer
        self.avatars = {}
        self.max_connections = 100

    def requestAvatar(self, avatarID, mind, *interfaces):
        if not pb.IPerspective in interfaces:
            raise NotImplementedError("No supported avatar interface.")
        else:
            if avatarID == 'admin':
                kind = 'admin'
            elif avatarID == 'Anonymous' and mind:
                kind = 'anonymous_monitor'
            elif mind:
                kind = 'monitor'
            else:
                kind = 'client'

            if (avatarID, kind) in self.avatars.keys():
                avatar = self.avatars[(avatarID, kind)]
            else:
                if kind == 'admin':
                    avatar = AdminPerspective(self.DSageServer, avatarID)
                elif kind == 'anonymous_monitor':
                    avatar = AnonymousMonitorPerspective(self.DSageServer,
                                                         avatarID)
                elif kind == 'monitor':
                    avatar = MonitorPerspective(self.DSageServer, avatarID)
                elif kind == 'client':
                    avatar = UserPerspective(self.DSageServer, avatarID)
                self.avatars[(avatarID, kind)] = avatar
        if avatar.connections >= self.max_connections:
            raise ValueError('Too many connections for user %s' % avatarID)

        avatar.attached(avatar, mind)
        log.msg('(%s, %s) connected' % (avatarID, kind))

        return pb.IPerspective, avatar, lambda a=avatar:a.detached(avatar,
                                                                   mind)