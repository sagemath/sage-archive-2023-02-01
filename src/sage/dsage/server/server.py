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
import datetime
import xml.dom.minidom
import cStringIO

from twisted.spread import pb
from twisted.python import log

from sage.dsage.errors.exceptions import BadTypeError
from sage.dsage.database.job import expand_job
from sage.dsage.misc.misc import timedelta_to_seconds


class DSageServer(pb.Root):
    """
    This class represents Distributed Sage server which does all the
    coordination of distributing jobs, creating new jobs and accepting job
    submissions.

    """

    def __init__(self, jobdb, monitordb, clientdb, log_level=0):
        """
        Initializes the Distributed Sage PB Server.

        Parameters:
        jobdb -- pass in the Database object
        log_level -- specifies the amount of logging to be done (default=0)

        """

        self.jobdb = jobdb
        self.monitordb = monitordb
        self.clientdb = clientdb
        self.LOG_LEVEL = log_level

    def register_client_factory(self, client_factory):
        self.client_factory = client_factory

    def unpickle(self, pickled_job):
        return cPickle.loads(zlib.decompress(pickled_job))

    def get_job(self, anonymous=False):
        """
        Returns a job to the client.

        This method returns the first job that has not been completed
        in our job database.

        """

        if anonymous:
            jdict = self.jobdb.get_job(anonymous=True)
        else:
            jdict = self.jobdb.get_job(anonymous=False)
        if jdict == None:
            if self.LOG_LEVEL > 3:
                log.msg('[DSage, get_job]' + ' Job db is empty.')
            return None
        else:
            job_id = jdict['job_id']
            if self.LOG_LEVEL > 3:
                log.msg('[DSage, get_job]' + ' Sending job %s' % job_id)
            jdict['status'] = 'processing'
            jdict['start_time'] = datetime.datetime.now()
            self.jobdb.store_job(jdict)

        return jdict

    def set_job_uuid(self, job_id, uuid):
        return self.jobdb.set_job_uuid(job_id, uuid)

    def set_busy(self, uuid, busy):
        return self.monitordb.set_busy(uuid, busy=busy)

    def get_job_by_id(self, job_id):
        """
        Returns a job by the job id.

        Parameters:
        id -- the job id

        """

        jdict = self.jobdb.get_job_by_id(job_id)
        uuid = jdict['monitor_id']
        jdict['worker_info'] = self.monitordb.get_monitor(uuid)

        return jdict

    def get_job_result_by_id(self, job_id):
        """Returns the job result.

        Parameters:
        id -- the job id (str)

        """

        jdict = self.jobdb.get_job_by_id(job_id)
        job = expand_job(jdict)

        return job.result

    def get_job_output_by_id(self, job_id):
        """Returns the job output.

        Parameters:
        id -- the job id (str)

        """

        job = self.jobdb.get_job_by_id(job_id)

        return job.output

    def sync_job(self, job_id):
        raise NotImplementedError

    def get_jobs_by_username(self, username, active=False):
        """
        Returns jobs created by username.

        Parameters:
        username -- the username (str)
        active -- when set to True, only return active jobs (bool)
        job_name -- the job name (optional)

        """

        jobs = self.jobdb.get_jobs_by_username(username, active)

        if self.LOG_LEVEL > 3:
            log.msg(jobs)

        return jobs

    def submit_job(self, jdict):
        """
        Submits a job to the job database.

        Parameters:
        jdict -- the internal dictionary of a Job object

        """

        if self.LOG_LEVEL > 3:
            log.msg('[DSage, submit_job] %s' % (jdict))
        if jdict['code'] is None:
            return False
        if jdict['name'] is None:
            jdict['name'] = 'Default'

        # New jobs should not have a job_id
        jdict['job_id'] = None

        jdict['update_time'] = datetime.datetime.now()

        job_id = self.jobdb.store_job(jdict)
        log.msg('Received job %s' % job_id)

        return job_id

    def get_all_jobs(self):
        """
        Returns a list of all jobs in the database.

        """
        return self.jobdb.get_all_jobs()

    def get_active_jobs(self):
        """
        Returns a list of active jobs"""

        return self.jobdb.get_active_jobs()

    def get_active_clients_list(self):
        """
        Returns a list of active clients.

        """

        raise NotImplementedError

    def get_killed_jobs_list(self):
        """
        Returns a list of killed job jdicts.

        """

        killed_jobs = self.jobdb.get_killed_jobs_list()

        return killed_jobs

    def get_next_job_id(self):
        """
        Returns the next job id.

        """

        if self.LOG_LEVEL > 0:
            log.msg('[DSage, get_next_job_id] Returning next job ID')

        return self.jobdb.get_next_job_id()

    def job_done(self, job_id, output, result, completed, cpu_time):
        """
        job_done is called by the workers check_output method.

        Parameters:
        job_id -- job id (string)
        output -- the stdout from the worker (string)
        result -- the result from the client (compressed pickle string)
                  result could be 'None'
        completed -- whether or not the job is completed (bool)

        """

        if self.LOG_LEVEL > 0:
            log.msg('[DSage, job_done] %s called back' % (job_id))
        if self.LOG_LEVEL > 3:
            log.msg('[DSage, job_done] output: %s ' % output)
            log.msg('[DSage, job_done] completed: %s ' % completed)

        jdict = self.get_job_by_id(job_id)
        output = str(output)
        if jdict['output'] is not None: # Append new output to existing output
            jdict['output'] += output
        else:
            jdict['output'] = output
        if completed:
            jdict['result'] = result
            jdict['cpu_time'] = cpu_time
            delta = timedelta_to_seconds(datetime.datetime.now() -
                                         jdict['start_time'])
            jdict['wall_time'] = delta
            jdict['status'] = 'completed'
            log.msg('%s completed!' % job_id)

        jdict['update_time'] = datetime.datetime.now()

        return self.jobdb.store_job(jdict)

    def job_failed(self, job_id, traceback):
        """
        job_failed is called when a remote job fails.

        Parameters:
        job_id -- the job id (str)

        """

        job = expand_job(self.jobdb.get_job_by_id(job_id))
        job.failures += 1
        job.output = traceback

        if job.failures > self.jobdb.job_failure_threshold:
            job.status = 'failed'
        else:
            job.status = 'new' # Put job back in the queue

        if self.LOG_LEVEL > 1:
            s = ['[DSage, job_failed] Job %s failed ' % (job_id),
                 '%s times. ' % (job.failures)]
            log.msg(''.join(s))
            if job.status == 'failed':
                msg = '%s failed, removing from queue.' % (job_id)
                log.msg(msg)

        job.update_time = datetime.datetime.now()

        return self.jobdb.store_job(job.reduce())

    def kill_job(self, job_id, reason):
        """
        Kills a job.

        Marks as job as killed and moves it to the killed jobs database.

        """

        if job_id == None:
            if self.LOG_LEVEL > 0:
                log.msg('[DSage, kill_job] Invalid job id')
            return None
        else:
            try:
                self.jobdb.set_killed(job_id, killed=True)
                if self.LOG_LEVEL > 0:
                    log.msg('Killed job %s' % (job_id))
            except Exception, msg:
                log.err(msg)
                log.msg('Failed to kill job %s' % job_id)

        return job_id

    def get_monitor_list(self):
        """
        Returns a list of workers as a 3 tuple.

        tuple[0] = broker object
        tuple[1] = ip
        tuple[2] = port

        """

        return self.monitordb.get_monitor_list()

    def get_client_list(self):
        """
        Returns a list of clients.

        """

        return self.clientdb.get_client_list()

    def get_cluster_speed(self):
        """
        Returns an approximation of the total CPU speed of the cluster.

        """
        raise NotImplementedError


    def get_worker_count(self):
        """
        Returns a list of busy and free workers.

        """

        count = {}
        free_workers = self.monitordb.get_worker_count(connected=True,
                                                       busy=False)
        working_workers = self.monitordb.get_worker_count(connected=True,
                                                          busy=True)

        count['free'] = free_workers
        count['working'] = working_workers

        return count

class DSageWorkerServer(DSageServer):
    """
    Exposes methods to workers.
    """

    def remote_get_job(self):
        return DSageServer.get_job(self)

    def remote_job_done(self, job_id, output, result, completed, worker_info):
        if not (isinstance(job_id, str) or isinstance(completed, bool)):
            log.msg('BadType in remote_job_done')
            raise BadTypeError()

        return DSageServer.job_done(self, job_id, output, result,
                             completed, worker_info)

    def remote_job_failed(self, job_id, traceback):
        if not isinstance(job_id, str):
            log.msg('BadType in remote_job_failed')
            raise BadTypeError()

        return DSageServer.job_failed(self, job_id, traceback)

    def remote_get_killed_jobs_list(self):
        return DSageServer.get_killed_jobs_list(self)