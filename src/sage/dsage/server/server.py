##############################################################################
#
#  DSAGE: Distributed SAGE
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
##############################################################################

import zlib
import cPickle

from twisted.spread import pb
from twisted.python import log

from sage.dsage.misc.hostinfo import HostInfo
import sage.dsage.server.worker_tracker as worker_tracker
from sage.dsage.server.hostinfo_tracker import hostinfo_list
from sage.dsage.errors.exceptions import BadTypeError

pb.setUnjellyableForClass(HostInfo, HostInfo)

class DSageServer(pb.Root):
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

class DSageWorkerServer(DSageServer):
    r"""
    Exposes methods to workers.
    """

    def remote_getJob(self):
        return DSageServer.getJob(self)

    def remote_jobDone(self, jobID, output, result, completed, worker_info):
        if not (isinstance(jobID, str) or isinstance(completed, bool)):
            raise BadTypeError()

        return DSageServer.jobDone(self, jobID, output, result,
                             completed, worker_info)

    def remote_jobFailed(self, jobID):
        if not isinstance(jobID, str):
            raise BadTypeError()

        return DSageServer.jobFailed(self, jobID)

    def remote_getKilledJobsList(self):
        return DSageServer.getKilledJobsList(self)

    def remote_submitHostInfo(self, hostinfo):
        if not isinstance(hostinfo, dict):
            raise BadTypeError()
        return DSageServer.submitHostInfo(self, hostinfo)

