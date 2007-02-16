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


import datetime
import zlib
import cPickle

from twisted.spread import pb
from zope.interface import implements
from twisted.cred import portal, credentials
from twisted.spread.interfaces import IJellyable
from twisted.spread.pb import IPerspective, AsReferenceable
from twisted.python import log
from twisted.internet import defer

from sage.dsage.database.job import Job
from sage.dsage.misc.hostinfo import HostInfo
import sage.dsage.server.client_tracker as client_tracker
import sage.dsage.server.worker_tracker as worker_tracker
from sage.dsage.server.hostinfo_tracker import hostinfo_list
from sage.dsage.errors.exceptions import BadJobError, BadTypeError

pb.setUnjellyableForClass(HostInfo, HostInfo)

class DSage(pb.Root):
    r"""
    This class represents Distributed Sage server which does all the
    coordination of distributing jobs, creating new jobs and accepting job
    submissions.

    """
    def __init__(self, jobdb, log_level=0):
        r"""
        Initializes the Distributed Sage PB Server.

        Parameters:
        jobdb -- pass in the Database object
        log_level -- specifies the amount of logging to be done (default=0)

        """

        self.jobdb = jobdb
        self.log_level = log_level

    def unpickle(self, pickled_job):
        return cPickle.loads(zlib.decompress(pickled_job))

    def getJob(self):
        r"""
        Returns a job to the client.

        This method returns the first job that has not been completed
        in our job database.

        """

        job = self.jobdb.getJob()

        if job == None:
            if self.log_level > 2:
                log.msg('[DSage]' + ' Job db is empty.')
            return None
        else:
            if self.log_level > 0:
                log.msg('[DSage, getJob]' + ' Returning Job (%s, %s) to client'
                % (job.id, job.name))
                job.status = 'processing'
            self.jobdb.storeJob(job)
        return job.pickle()

    def getJobByID(self, jobID):
        r"""
        Returns a job by the job id.

        Parameters:
        id -- the job id

        """

        job = self.jobdb.getJobByID(jobID)
        return job.pickle()

    def getJobResultByID(self, jobID):
        """Returns the job result.

        Parameters:
        id -- the job id (str)

        """

        job = self.jobdb.getJobByID(jobID)
        return job.result

    def getJobOutputByID(self, jobID):
        """Returns the job output.

        Parameters:
        id -- the job id (str)

        """

        job = self.jobdb.getJobByID(jobID)

        return job.output

    def syncJob(self, jobID):
        job = self.jobdb.getJobByID(jobID)
        # new_job = copy.deepcopy(job)
        # print new_job
        # # Set file, data to 'Omitted' so we don't need to transfer it
        # new_job.file = 'Omitted...'
        # new_job.data = 'Omitted...'

        return job.pickle()

    def getJobsByAuthor(self, author, is_active, job_name):
        r"""
        Returns jobs created by author.

        Parameters:
        author -- the author name (str)
        is_active -- when set to True, only return active jobs (bool)
        job_name -- the job name (optional)

        """

        jobs = self.jobdb.getJobsByAuthor(author, is_active, job_name)

        if self.log_level > 3:
            log.msg(jobs)
        return jobs

    def submitJob(self, pickled_job):
        r"""
        Submits a job to the job database.

        Parameters:
        pickled_job -- a pickled_Job object

        """

        try:
            # job = cPickle.loads(zlib.decompress(pickled_job))
            # print job
            if self.log_level > 3:
                log.msg('[DSage, submitJob] Trying to unpickle job.')
            job = self.unpickle(pickled_job)
            # print job
        except:
            return False
        if job.file is None:
            return False
        if job.name is None:
            job.name = 'defaultjob'
        if job.id is None:
            job.id = self.jobdb.getJobID()

        if self.log_level > 0:
            log.msg('[DSage, submitJob] Job (%s %s) submitted' % (job.id,
                                                                 job.name))
        return self.jobdb.storeJob(job)

    def getJobsList(self):
        r"""
        Returns an ordered list of jobs in the database.

        """
        return self.jobdb.getJobsList()

    def getActiveJobs(self):
        r"""
        Returns a list of active jobs"""

        return self.jobdb.getActiveJobs()

    def getActiveClientsList(self):
        r"""
        Returns a list of active clients.

        """

        raise NotImplementedError

    def getKilledJobsList(self):
        r"""
        Returns a list of killed jobs.
        """

        killed_jobs = self.jobdb.getKilledJobsList()
        return [job.pickle() for job in killed_jobs]

    def getNextJobID(self):
        r"""
        Returns the next job id.

        """

        if self.log_level > 0:
            log.msg('[DSage, getNextJobID] Returning next job ID')

        return self.jobdb.getJobID()

    def jobDone(self, jobID, output, result, completed, worker_info):
        r"""
        jobDone is called by the workers checkForJobOutput method.

        Parameters:
        jobID -- job id (str)
        output -- the stdout from the worker (string)
        result -- the result from the client (compressed pickle string)
                  result could be 'None'
        completed -- whether or not the job is completed (bool)
        worker_info -- ''.join(os.uname())

        """

        if self.log_level > 0:
            log.msg('[DSage, jobDone] Job %s called back' % (jobID))
        if self.log_level > 2:
            log.msg('[DSage, jobDone] Output: %s ' % output)
            log.msg('[DSage, jobDone] Result: Some binary data...')
            log.msg('[DSage, jobDone] completed: %s ' % completed)
            log.msg('[DSage, jobDone] worker_info: %s ' % str(worker_info))

        job = self.unpickle(self.getJobByID(jobID))

        if self.log_level > 3:
            log.msg('[DSage, jobDone] result type' , type(result))

        output = str(output)
        if job.output is not None: # Append new output to existing output
            job.output = job.output + output
        else:
            job.output = output
        if completed:
            job.result = result
            job.status = 'completed'
            job.worker_info = worker_info

        return self.jobdb.storeJob(job)

    def jobFailed(self, jobID):
        r"""
        jobFailed is called when a remote job fails.

        Parameters:
        jobID -- the job id (str)

        """

        job = self.jobdb.getJobByID(jobID)
        job.failures += 1
        job.status = 'new' # Put job back in the queue

        if self.log_level > 1:
            log.msg('[DSage, jobFailed] Job %s failed' % jobID)

        self.jobdb.storeJob(job)

    def killJob(self, jobID, reason):
        r"""
        Kills a job.

        Marks as job as killed and moves it to the killed jobs database.

        """

        job = self.unpickle(self.getJobByID(jobID))
        if job == None:
            if self.log_level > 0:
                log.msg('[DSage, killJob] No such job id')
            return None
        else:
            job.killed = True
            self.jobdb.storeJob(job)
            if self.log_level > 0:
                log.msg('Job %s was killed because %s ' % (jobID, reason))

        return jobID

    def getWorkersList(self):
        r"""
        Returns a list of workers as a 3 tuple.

        tuple[0] = broker object
        tuple[1] = ip
        tuple[2] = port

        """

        return worker_tracker.worker_list

    def getClusterSpeed(self):
        r"""
        Returns an approximation of the total CPU speed of the cluster.

        """

        cluster_speed = 0
        if self.log_level > 3:
            log.msg(hostinfo_list)
            log.msg(len(hostinfo_list))
        for h in hostinfo_list:
            speed_multiplier = int(h['cpus'])
            for k,v in h.iteritems():
                if k == 'cpu_speed':
                    cluster_speed += float(v) * speed_multiplier

        return cluster_speed

    def submitHostInfo(self, h):
        r"""
        Takes a dict of workers machine specs.

        """

        if self.log_level > 0:
            log.msg(h)
        if len(hostinfo_list) == 0:
            hostinfo_list.append(h)
        else:
            for h in hostinfo_list:
                if h['uuid'] not in h.values():
                    hostinfo_list.append(h)

class DSageForWorkers(DSage):
    r"""
    Exposes methods to workers.
    """

    def remote_getJob(self):
        return DSage.getJob(self)

    def remote_jobDone(self, jobID, output, result, completed, worker_info):
        if not (isinstance(jobID, str) or isinstance(completed, bool)):
            raise BadTypeError()

        return DSage.jobDone(self, jobID, output, result,
                             completed, worker_info)

    def remote_jobFailed(self, jobID):
        if not isinstance(jobID, str):
            raise BadTypeError()

        return DSage.jobFailed(self, jobID)

    def remote_getKilledJobsList(self):
        return DSage.getKilledJobsList(self)

    def remote_submitHostInfo(self, hostinfo):
        if not isinstance(hostinfo, dict):
            raise BadTypeError()
        return DSage.submitHostInfo(self, hostinfo)

class WorkerPBServerFactory(pb.PBServerFactory):
    def __init__(self, root):
        pb.PBServerFactory.__init__(self, root)

    def clientConnectionMade(self, broker):
        """Keeps a 3-tuple of connected workers.
        tuple[0] - the broker
        tuple[1] - the worker ip
        tuple[2] - the worker port

        """
        # self.clients.append((broker.transport.protocol,
        #                      broker.transport.getPeer().host,
        #                      broker.transport.getPeer().port))

        worker_tracker.add((broker,
                            broker.transport.getPeer().host,
                            broker.transport.getPeer().port))


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

    def __init__(self, DSage):
        self.DSage = DSage

    def perspective_getJob(self):
        return self.DSage.getJob()

    def perspective_getNextJobID(self):
        job_id = self.DSage.getNextJobID()

        return job_id

    def perspective_getJobByID(self, jobID):
        if not isinstance(jobID, str):
            raise BadTypeError()

        job = self.DSage.getJobByID(jobID)

        return job

    def perspective_getJobsByAuthor(self, author, job_name, is_active):
        if not (isinstance(author, str) or
                isinstance(is_active, bool) or
                isinstance(job_name, str)):
            raise BadTypeError()

        jobs = self.DSage.getJobsByAuthor(author, job_name, is_active)

        return jobs

    def perspective_getJobResultByID(self, jobID):
        if not isinstance(jobID, str):
            raise BadTypeError()

        return self.DSage.getJobResultByID(jobID)

    def perspective_getJobOutputByID(self, jobID):
        if not isinstance(jobID, str):
            raise BadTypeError()

        return self.DSage.getJobOutputByID(jobID)

    def perspective_syncJob(self, jobID):
        if not isinstance(jobID, str):
            return None

        return self.DSage.syncJob(jobID)

    def perspective_submitJob(self, pickled_job):
        return self.DSage.submitJob(pickled_job)

    def perspective_jobDone(self, jobID, output, result,
                            completed, worker_info):
        if not (isinstance(jobID, str) or isinstance(completed, bool)):
            raise BadTypeError()

        return self.DSage.jobDone(jobID, output, result,
                                  completed, worker_info)

    def perspective_jobFailed(self, jobID):
        if not isinstance(jobID, str):
            raise BadTypeError()

        return self.DSage.jobFailed(jobID)

    def perspective_killJob(self, jobID, reason=None):
        if not isinstance(jobID, str):
            raise BadTypeError()

        return self.DSage.killJob(jobID, reason)

    def perspective_getClusterSpeed(self):
        return self.DSage.getClusterSpeed()

class AdminPerspective(UserPerspective):
    r"""
    Defines the perspective of the admin.

    """

    def __init__(self, DSage):
        UserPerspective.__init__(self, DSage)
        log.msg('[Realm]' + " admin connected")

    def perspective_getClientList(self):
        return client_tracker.client_list

    def perspective_getWorkerList(self):
        r"""
        Returns a list of 2-tuples of connected workers.

        """

        return [x[1:] for x in self.DSage.getWorkersList()]


class Realm(object):
    implements(portal.IRealm)

    def __init__(self, DSage):
        self.DSage = DSage

    def requestAvatar(self, avatarId, mind, *interfaces):
        if not pb.IPerspective in interfaces:
            raise NotImplementedError, "No supported avatar interface."
        else:
            if avatarId == 'admin':
                avatar = AdminPerspective(self.DSage)
            else:
                avatar = UserPerspective(self.DSage)
        avatar.attached(avatar, mind)
        return pb.IPerspective, avatar, lambda a=avatar:a.detached(avatar,
                                                                   mind)
