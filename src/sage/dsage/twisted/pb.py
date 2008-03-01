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
#
############################################################################

from twisted.spread import pb
from twisted.spread import banana
banana.SIZE_LIMIT = 100*1024*1024 # 100 MegaBytes

from zope.interface import implements

from twisted.internet import reactor
from twisted.cred import portal, credentials
from twisted.cred.credentials import ISSHPrivateKey
from twisted.cred.credentials import IAnonymous
from twisted.cred.credentials import Anonymous
from twisted.internet.protocol import ReconnectingClientFactory
from twisted.spread.interfaces import IJellyable
from twisted.spread.pb import (IPerspective, AsReferenceable, PBClientFactory)
from twisted.python import log

from sage.dsage.errors.exceptions import BadTypeError, BadJobError
from sage.dsage.misc.misc import gen_uuid, check_uuid

class ClientFactory(PBClientFactory, ReconnectingClientFactory):
    """
    Custom implementation of the ClientFactory that supports logging in
    with public key as well as unauthenticated credentials.

    """

    def __init__(self, cb, args, kwargs):
        PBClientFactory.__init__(self)
        self._observer = (cb, args, kwargs)

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
        elif IAnonymous.providedBy(creds):
            d = self.getRootObject()
            d.addCallback(self._cbAnonymousLogin, mind)
            return d
        else:
            raise TypeError('Invalid credentials.')

    def _cbSendUsername(self, root, username, algorithm, blob, sig_data,
                        signature, mind):
        d = root.callRemote("login", username, algorithm, blob, sig_data,
                                signature, mind)
        return d

    def _cbAnonymousLogin(self, root, mind):
        d = root.callRemote("login_authenticate", mind)

        return d

    def startConnecting(self, server, port, cred=None):
        if cred:
            reactor.connectTLS(self.server, self.port, self.factory,
                               cred)
        else:
            reactor.connectTCP(self.server, self.port, self.factory)

    def clientConnectionMade(self, broker):
        PBClientFactory.clientConnectionMade(self, broker)
        cb, args, kwargs = self._observer
        cb(self._root, *args, **kwargs)

    def clientConnectionLost(self, connector, reason, reconnecting=0):
        PBClientFactory.clientConnectionLost(self, connector, reason,
                                             reconnecting=1)
        ReconnectingClientFactory.clientConnectionLost(self, connector,
                                                       reason)

    def clientConnectionFailed(self, connector, reason):
        PBClientFactory.clientConnectionFailed(self, connector, reason)
        ReconnectingClientFactory.clientConnectionFailed(self, connector,
                                                         reason)

class _SSHKeyPortalRoot(pb._PortalRoot):
    def rootObject(self, broker):
        return _SSHKeyPortalWrapper(self.portal, broker)

class _SSHKeyPortalWrapper(pb._PortalWrapper):
    def remote_login(self, username, algorithm, blob, data, signature, mind):
        pubkey_cred = credentials.SSHPrivateKey(username,
                                                algorithm,
                                                blob,
                                                data,
                                                signature)

        d = self.portal.login(pubkey_cred, mind, IPerspective)
        d.addCallback(self._loggedIn)

        return d

    def remote_login_authenticate(self, mind):
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

    def __init__(self, dsage_server, avatarID):
        self.dsage_server = dsage_server
        self.avatarID = avatarID
        self.connections = 0
        self.workerdb = self.dsage_server.workerdb
        self.clientdb = self.dsage_server.clientdb
        self.kind = ''
        self.host_info = []

    def __repr__(self):
        return "<%s:'%s'>" % (self.kind, self.avatarID)

    def perspectiveMessageReceived(self, broker, message, args, kw):
        self.broker = broker

        return pb.Avatar.perspectiveMessageReceived(self, broker,
                                                    message, args, kw)


    def attached(self, avatar, mind):
        self.connections += 1
        log.msg('%s connected' % avatar)


    def detached(self, avatar, mind):
        self.connections -= 1
        log.msg('%s disconnected' % avatar)



class AnonymousWorker(DefaultPerspective):
    """
    Defines the perspective of an authenticate worker.

    """

    def __init__(self, dsage_server, avatarID):
        DefaultPerspective.__init__(self, dsage_server, avatarID)

    def attached(self, avatar, mind):
        DefaultPerspective.attached(self, avatar, mind)
        self.monitor = mind[0]
        self.host_info = mind[1]
        self.host_info['ip'] = mind[0].broker.transport.getPeer().host
        self.host_info['port'] = mind[0].broker.transport.getPeer().port
        if check_uuid(self.host_info['uuid']):
            uuid = self.host_info['uuid']
            if self.workerdb.get_worker(uuid) is None:
                self.workerdb.add_worker(self.host_info)
            else:
                self.workerdb.update_worker(self.host_info)
        else:
            uuid = gen_uuid()
            if isinstance(self.monitor, pb.RemoteReference):
                d = self.monitor.callRemote('set_uuid', uuid)
                self.host_info['uuid'] = uuid
                self.workerdb.add_worker(self.host_info)
        self.uuid = uuid
        self.workerdb.set_connected(uuid, connected=True)
        self.workerdb.set_authenticated(uuid, False)
        self.dsage_server.workers[uuid] = mind[0]

        return uuid

    def detached(self, avatar, mind):
        DefaultPerspective.detached(self, avatar, mind)
        self.workerdb.set_connected(self.uuid, connected=False)
        del self.dsage_server.workers[self.uuid]

    def perspective_get_job(self):
        """
        Returns jobs only marked as doable by authenticate workers.

        """

        uuid = self.host_info['uuid']
        jdict = self.dsage_server.get_job()
        if jdict is not None:
            self.dsage_server.set_job_uuid(jdict['job_id'], uuid)
            self.dsage_server.set_busy(uuid, True)
        else:
            self.dsage_server.set_busy(uuid, False)

        return jdict

    def perspective_get_killed_jobs_list(self):
        return self.dsage_server.get_killed_jobs_list()

    def perspective_job_failed(self, job_id, traceback):
        if not isinstance(job_id, str):
            log.msg('Bad job_id %s' % (job_id))
            raise BadTypeError()

        uuid = self.host_info['uuid']
        self.dsage_server.set_busy(uuid, False)

        return self.dsage_server.job_failed(job_id, traceback)

    def perspective_job_done(self, job_id, output, result, completed,
                             cpu_time):
        if not (isinstance(job_id, str) or isinstance(completed, bool)):
            log.msg('Bad job_id passed to perspective_job_done')
            log.msg('job_id: %s' % (job_id))
            log.msg('output: %s' % (output))
            log.msg('completed: %s' % (completed))
            raise BadTypeError()
        if completed:
            uuid = self.host_info['uuid']
            self.dsage_server.set_busy(uuid, False)

        return self.dsage_server.job_done(job_id, output, result, completed,
                                         cpu_time)

class Worker(AnonymousWorker):
    """
    Defines the perspective of an authenticated worker to the server.

    """
    def __init__(self, dsage_server, avatarID):
        DefaultPerspective.__init__(self, dsage_server, avatarID)

    def attached(self, avatar, mind):
        uuid = AnonymousWorker.attached(self, avatar, mind)
        self.workerdb.set_authenticated(uuid, True)
        return uuid

    def perspective_get_job(self):
        """
        Returns jobs to authenticated workers.

        """

        try:
            uuid = self.host_info['uuid']
        except Exception, msg:
            raise ValueError("Could not match a uuid to the monitor.")

        jdict = self.dsage_server.get_job()
        if jdict is not None:
            self.dsage_server.set_job_uuid(jdict['job_id'], uuid)
            self.dsage_server.set_busy(uuid, True)

        return jdict

class Client(DefaultPerspective):
    """
    Defines the perspective of a regular user to the server.

    """

    def __init__(self, dsage_server, avatarID):
        DefaultPerspective.__init__(self, dsage_server, avatarID)

    def attached(self, avatar, mind):
        DefaultPerspective.attached(self, avatar, mind)
        self.clientdb.set_connected(self.avatarID, connected=True)
        self.clientdb.update_login_time(self.avatarID)
        self.dsage_server.clients.append(self)

    def detached(self, avatar, mind):
        DefaultPerspective.detached(self, avatar, mind)
        self.clientdb.set_connected(self.avatarID, connected=False)
        self.dsage_server.clients.remove(self)

    def perspective_get_job_by_id(self, job_id):
        if not isinstance(job_id, str):
            log.msg('Bad job_id [%s] passed to get_job_by_id' % (job_id))
            raise BadTypeError()
        job = self.dsage_server.get_job_by_id(job_id)

        return job

    def perspective_get_jobs_by_username(self, username, status):
        if not (isinstance(username, str)):
            log.msg('Bad username [%s] passed to ' +
                    'perspective_get_jobs_by_username' % (username))
            raise BadTypeError()

        jobs = self.dsage_server.get_jobs_by_username(username, status)

        return jobs

    def perspective_get_job_result_by_id(self, job_id):
        if not isinstance(job_id, str):
            log.msg('Bad job_id [%s] passed to' +
                    'perspective_get_job_result_by_id' % (job_id))
            raise BadTypeError()

        return self.dsage_server.get_job_result_by_id(job_id)

    def perspective_get_job_output_by_id(self, job_id):
        if not isinstance(job_id, str):
            log.msg('Bad job_id [%s] passed to ' +
                    'get_job_output_by_id' % (job_id))
            raise BadTypeError()

        return self.dsage_server.get_job_output_by_id(job_id)

    def perspective_sync_job(self, job_id):
        if not isinstance(job_id, str):
            return None

        return self.dsage_server.sync_job(job_id)

    def perspective_submit_job(self, jdict):
        if jdict is None:
            raise BadJobError()

        return self.dsage_server.submit_job(jdict)

    def perspective_web_server_url(self):
        return "http://localhost:%s" % self.dsage_server.web_port

    def perspective_kill_job(self, job_id):
        if not isinstance(job_id, str):
            log.msg('Bad job_id [%s] passed to perspective_kill_job' % job_id)
            raise BadTypeError()

        return self.dsage_server.kill_job(job_id)

    def perspective_get_cluster_speed(self):
        return self.dsage_server.get_cluster_speed()

    def perspective_get_worker_list(self):
        # return [x[1] for x in self.dsage_server.get_worker_list()]
        return self.dsage_server.get_worker_list()

    def perspective_get_client_list(self):
        return self.dsage_server.get_client_list()

    def perspective_get_worker_count(self):
        return self.dsage_server.get_worker_count()

    def perspective_get_killed_jobs_list(self):
        return self.dsage_server.get_killed_jobs_list()

    def perspective_read_log(self, n, kind):
        return self.dsage_server.read_log(n, kind)


class Admin(Client, Worker):
    """
    Defines the perspective of the admin.

    """

    def __init__(self, dsage_server, avatarID):
        Client.__init__(self, dsage_server, avatarID)

class Tester(Client, Worker):
    def __init__(self, dsage_server, avatarID):
        DefaultPerspective.__init__(self, dsage_server, avatarID)

    def attached(self, avatar, mind):
        try:
            Worker.attached(self, avatar, mind)
        except:
            Client.attached(self, avatar, mind)

    def detached(self, avatar, mind):
        try:
            Worker.detached(self, avatar, mind)
        except:
            Client.attached(self, avatar, mind)

class Realm(object):
    implements(portal.IRealm)

    def __init__(self, dsage_server):
        self.dsage_server = dsage_server
        self.max_connections = 1000
        self.client_avatars = {}

    def requestAvatar(self, avatarID, mind, *interfaces):
        if not pb.IPerspective in interfaces:
            raise NotImplementedError("No supported avatar interface.")
        else:
            if self.dsage_server._testing:
                kind = 'tester'
                avatarID = 'tester'
                avatar = Tester(self.dsage_server, avatarID)
            else:
                if avatarID == 'admin':
                    kind = 'admin'
                elif avatarID == 'Anonymous':
                    kind = 'unauthenticated_worker'
                elif mind:
                    kind = 'worker'
                else:
                    kind = 'client'
                if kind == 'admin':
                    avatar = Admin(self.dsage_server, avatarID)
                elif kind == 'client':
                    if avatarID in self.client_avatars:
                        avatar = self.client_avatars[avatarID]
                    else:
                        avatar = Client(self.dsage_server, avatarID)
                        self.client_avatars[avatarID] = avatar
                elif kind == 'unauthenticated_worker':
                    avatar = AnonymousWorker(self.dsage_server, avatarID)
                elif kind == 'worker':
                    avatar = Worker(self.dsage_server, avatarID)
        avatar.kind = kind
        avatar.attached(avatar, mind)
        self.max_connections += 1
        if avatar.connections >= self.max_connections:
            raise ValueError('Too many connections.')

        return pb.IPerspective, avatar, lambda a = avatar:a.detached(avatar,
                                                                     mind)
