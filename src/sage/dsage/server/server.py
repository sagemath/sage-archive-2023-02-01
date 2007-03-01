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
import sage.dsage.server.client_tracker as client_tracker
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

    def get_job(self):
        r"""
        Returns a job to the client.

        This method returns the first job that has not been completed
        in our job database.

        """

        job = self.jobdb.get_job()

        if job == None:
            if self.log_level > 2:
                log.msg('[DSage, get_job]' + ' Job db is empty.')
            return None
        else:
            if self.log_level > 3:
                log.msg('[DSage, get_job]' + ' Returning Job %s to client'
                        % (job.id))
                job.status = 'processing'
            self.jobdb.store_job(job)
        return job.pickle()

    def get_job_by_id(self, jobID):
        r"""
        Returns a job by the job id.

        Parameters:
        id -- the job id

        """

        job = self.jobdb.get_job_by_id(jobID)
        return job.pickle()

    def get_job_result_by_id(self, jobID):
        """Returns the job result.

        Parameters:
        id -- the job id (str)

        """

        job = self.jobdb.get_job_by_id(jobID)
        return job.result

    def get_job_output_by_id(self, jobID):
        """Returns the job output.

        Parameters:
        id -- the job id (str)

        """

        job = self.jobdb.get_job_by_id(jobID)

        return job.output

    def sync_job(self, jobID):
        job = self.jobdb.get_job_by_id(jobID)
        # new_job = copy.deepcopy(job)
        # print new_job
        # # Set file, data to 'Omitted' so we don't need to transfer it
        # new_job.file = 'Omitted...'
        # new_job.data = 'Omitted...'

        return job.pickle()

    def get_jobs_by_author(self, author, is_active, job_name):
        r"""
        Returns jobs created by author.

        Parameters:
        author -- the author name (str)
        is_active -- when set to True, only return active jobs (bool)
        job_name -- the job name (optional)

        """

        jobs = self.jobdb.get_jobs_by_author(author, is_active, job_name)

        if self.log_level > 3:
            log.msg(jobs)
        return jobs

    def submit_job(self, pickled_job):
        r"""
        Submits a job to the job database.

        Parameters:
        pickled_job -- a pickled_Job object

        """

        try:
            if self.log_level > 3:
                log.msg('[DSage, submit_job] Trying to unpickle job')
            job = self.unpickle(pickled_job)
        except:
            return False
        if job.file is None:
            return False
        if job.name is None:
            job.name = 'defaultjob'
        if job.id is None:
            job.id = self.jobdb.get_next_job_id()

        if self.log_level > 0:
            log.msg('[DSage, submit_job] Job (%s %s) submitted' % (job.id,
                                                                 job.name))
        return self.jobdb.store_job(job)

    def get_jobs_list(self):
        r"""
        Returns an ordered list of jobs in the database.

        """
        return self.jobdb.get_jobs_list()

    def get_active_jobs(self):
        r"""
        Returns a list of active jobs"""

        return self.jobdb.get_active_jobs()

    def get_active_clients_list(self):
        r"""
        Returns a list of active clients.

        """

        raise NotImplementedError

    def get_killed_jobs_list(self):
        r"""
        Returns a list of killed jobs.
        """

        killed_jobs = self.jobdb.get_killed_jobs_list()
        return [job.pickle() for job in killed_jobs]

    def get_next_job_id(self):
        r"""
        Returns the next job id.

        """

        if self.log_level > 0:
            log.msg('[DSage, get_next_job_id] Returning next job ID')

        return self.jobdb.get_next_job_id()

    def job_done(self, jobID, output, result, completed, worker_info):
        r"""
        job_done is called by the workers checkForJobOutput method.

        Parameters:
        jobID -- job id (str)
        output -- the stdout from the worker (string)
        result -- the result from the client (compressed pickle string)
                  result could be 'None'
        completed -- whether or not the job is completed (bool)
        worker_info -- ''.join(os.uname())

        """

        if self.log_level > 0:
            log.msg('[DSage, job_done] Job %s called back' % (jobID))
        if self.log_level > 2:
            log.msg('[DSage, job_done] Output: %s ' % output)
            log.msg('[DSage, job_done] Result: Some binary data...')
            log.msg('[DSage, job_done] completed: %s ' % completed)
            log.msg('[DSage, job_done] worker_info: %s ' % str(worker_info))

        job = self.unpickle(self.get_job_by_id(jobID))

        if self.log_level > 3:
            log.msg('[DSage, job_done] result type' , type(result))

        output = str(output)
        if job.output is not None: # Append new output to existing output
            job.output = job.output + output
        else:
            job.output = output
        if completed:
            job.result = result
            job.status = 'completed'
            job.worker_info = worker_info

        return self.jobdb.store_job(job)

    def job_failed(self, jobID):
        r"""
        job_failed is called when a remote job fails.

        Parameters:
        jobID -- the job id (str)

        """

        job = self.jobdb.get_job_by_id(jobID)
        job.failures += 1

        if job.failures > self.jobdb.JOB_FAILURE_THRESHOLD:
            job.status = 'failed'
        else:
            job.status = 'new' # Put job back in the queue

        if self.log_level > 1:
            s = ['[DSage, job_failed] Job %s failed ' % (jobID),
                 '%s times. ' % (job.failures)]
            log.msg(''.join(s))
        self.jobdb.store_job(job)

    def kill_job(self, jobID, reason):
        r"""
        Kills a job.

        Marks as job as killed and moves it to the killed jobs database.

        """

        job = self.unpickle(self.get_job_by_id(jobID))
        if job == None:
            if self.log_level > 0:
                log.msg('[DSage, kill_job] No such job id')
            return None
        else:
            job.killed = True
            self.jobdb.store_job(job)
            if self.log_level > 0:
                log.msg('Job %s was killed because %s ' % (jobID, reason))

        return jobID

    def get_worker_list(self):
        r"""
        Returns a list of workers as a 3 tuple.

        tuple[0] = broker object
        tuple[1] = ip
        tuple[2] = port

        """

        return worker_tracker.worker_list

    def get_client_list(self):
        r"""
        Returns a list of clients.

        """

        return client_tracker.client_list

    def get_cluster_speed(self):
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

    def submit_host_info(self, h):
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

    def remote_get_job(self):
        return DSageServer.get_job(self)

    def remote_job_done(self, jobID, output, result, completed, worker_info):
        if not (isinstance(jobID, str) or isinstance(completed, bool)):
            log.msg('BadType in remote_job_done')
            raise BadTypeError()

        return DSageServer.job_done(self, jobID, output, result,
                             completed, worker_info)

    def remote_job_failed(self, jobID):
        if not isinstance(jobID, str):
            log.msg('BadType in remote_job_failed')
            raise BadTypeError()

        return DSageServer.job_failed(self, jobID)

    def remote_get_killed_jobs_list(self):
        return DSageServer.get_killed_jobs_list(self)

    def remote_submit_host_info(self, hostinfo):
        if not isinstance(hostinfo, dict):
            log.msg('BadType in remote_submit_host_info')
            raise BadTypeError()
        return DSageServer.submit_host_info(self, hostinfo)

