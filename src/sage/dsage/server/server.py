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
import datetime
import os

from twisted.spread import pb
from twisted.python import log
from subprocess import Popen

from sage.dsage.database.job import expand_job
from sage.dsage.misc.misc import timedelta_to_seconds
from sage.dsage.misc.constants import SERVER_LOG, WORKER_LOG

class DSageServer(pb.Root):
    """
    Distributed Sage server which does all the coordination of distributing
    jobs, creating new jobs and accepting job submissions.

    """

    def __init__(self, jobdb, workerdb, clientdb, log_level=0):
        """
        Initializes the Distributed Sage PB Server.

        :type jobdb: sage.dsage.database.jobdb.JobDatabaseSQLite
        :param jobdb: a instance of the job database

        :type workerdb: sage.dsage.database.workerdb.WorkerDatabase
        :param workerdb: instance of the monitor database

        :type log_level: integer
        :param log_level: level of logging verbosity, higher is more verbose

        """

        self.jobdb = jobdb
        self.workerdb = workerdb
        self.clientdb = clientdb
        self.log_level = log_level
        self.clients = []
        self.workers = {}

        # Setting this to true results in NO authentication being made.
        self._testing = False

    def get_job(self):
        """
        Returns a job to the client.

        This method returns the first job that has not been completed
        in our job database.

        :type authenticated: boolean
        :param authenticated: whether or not the requester is authenticated

        """

        job = self.jobdb.get_job()
        if job == None:
            if self.log_level > 3:
                log.msg('[dsage, get_job]' + ' Job db is empty.')
            return None
        else:
            job_id = job.job_id
            if self.log_level > 3:
                log.msg('[dsage, get_job]' + ' Returning job %s' % job_id)
            job.status = 'processing'
            job.start_time = datetime.datetime.now()
            self.jobdb.update_job(job)

        return job._reduce()

    def set_job_uuid(self, job_id, uuid):
        """
        Sets the job's universal unique identifer, which identifies the worker
        that processed the job.

        :type job_id: string
        :param job_id: unique job identifier

        :type uuid: string
        :param uuid: universial unique identifier for the worker

        """

        return self.jobdb.set_job_uuid(job_id, uuid)

    def set_busy(self, uuid, busy):
        """
        Sets whether or not a particular worker is busy.

        :type uuid: string
        :param uuid: universial unique identifier

        :type busy: boolean
        :param busy: Whether or not the worker is busy

        """

        return self.workerdb.set_busy(uuid, busy=busy)

    def get_job_by_id(self, job_id):
        """
        Returns a job by the job id.

        :type job_id: string
        :param job_id: unique job identifier

        """

        jdict = self.jobdb.get_job_by_id(job_id)._reduce()

        return jdict

    def get_job_result_by_id(self, job_id):
        """Returns the job result.

        :type job_id: string
        :param job_id: unique job identifier

        """

        job = self.jobdb.get_job_by_id(job_id)

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

    def get_jobs_by_username(self, username, status):
        """
        Returns jobs created by username.

        Parameters:
        username -- the username (str)

        """

        jobs = self.jobdb.get_jobs_by_username(username, status)

        return [job._reduce() for job in jobs]

    def submit_job(self, jdict):
        """
        Submits a job to the job database.

        Parameters:
        jdict -- the internal dictionary of a Job object

        """

        if self.log_level > 3:
            log.msg('[submit_job] %s' % (jdict))
        if jdict['code'] is None:
            return False
        if jdict['name'] is None:
            jdict['name'] = 'Default'

        jdict['update_time'] = datetime.datetime.now()

        job_id = self.jobdb.store_jdict(jdict)
        log.msg('Received job %s' % job_id)

        self.push_job()

        return job_id

    def push_job(self):
        for worker in self.workers.values():
            worker.callRemote('get_job')

    def get_all_jobs(self):
        """
        Returns a list of all jobs in the database.

        """

        return [job._reduce() for job in self.jobdb.get_all_jobs()]

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

        return [job._reduce() for job in killed_jobs]

    def get_next_job_id(self):
        """
        Returns the next job id.

        """

        if self.log_level > 0:
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

        if self.log_level > 0:
            log.msg('[DSage, job_done] %s called back' % (job_id))
        if self.log_level > 3:
            log.msg('[DSage, job_done] output: %s ' % output)
            log.msg('[DSage, job_done] completed: %s ' % completed)

        output = str(output)
        job = self.jobdb.get_job_by_id(job_id)
        job.output += output
        job.wall_time = datetime.datetime.now() - job.start_time
        job.update_time = datetime.datetime.now()
        if completed:
            job.result = result
            job.cpu_time = cpu_time
            job.status = 'completed'
            job.finish_time = datetime.datetime.now()
        self.jobdb.sess.save_or_update(job)
        self.jobdb.sess.commit()

        return job_id

    def job_failed(self, job_id, traceback):
        """
        job_failed is called when a remote job fails.

        Parameters:
        job_id -- the job id (str)

        """

        job = self.jobdb.get_job_by_id(job_id)
        job.failures += 1
        job.output = traceback

        if job.failures > self.jobdb.failure_threshold:
            job.status = 'failed'
        else:
            job.status = 'new' # Put job back in the queue

        if self.log_level > 1:
            s = ['[DSage, job_failed] Job %s failed ' % (job_id),
                 '%s times. ' % (job.failures)]
            log.msg(''.join(s))
            if job.status == 'failed':
                msg = '%s failed, removing from queue.' % (job_id)
                log.msg(msg)

        job.update_time = datetime.datetime.now()

        return self.jobdb.store_jdict(job._reduce())

    def kill_job(self, job_id):
        """
        Kills a job.

        Marks as job as killed and moves it to the killed jobs database.

        """

        try:
            job = self.jobdb.set_killed(job_id, killed=True)
            if self.log_level > 0:
                log.msg('Killed job %s' % (job_id))
        except Exception, msg:
            log.err(msg)
            log.msg('Failed to kill job %s' % job_id)
            return None

        try:
            self.workers[job.uuid].callRemote('kill_job', job_id)
        except KeyError:
            pass

        return job_id

    def get_worker_list(self):
        """
        Returns a list of workers as a 3 tuple.

        tuple[0] = broker object
        tuple[1] = ip
        tuple[2] = port

        """

        return self.workerdb.get_worker_list()

    def get_client_list(self):
        """
        Returns a list of clients.

        """

        return [c.username for c in self.clientdb.get_client_list()]

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
        free_workers = self.workerdb.get_worker_count(connected=True,
                                                       busy=False)
        working_workers = self.workerdb.get_worker_count(connected=True,
                                                          busy=True)

        count['free'] = free_workers
        count['working'] = working_workers

        return count

    def upgrade_workers(self):
        """
        Upgrades the connected workers to the latest SAGE version.

        """

        raise NotImplementedError

    def read_log(self, n, kind):
        """
        Returns the last n lines of the server log.
        Defaults to returning the last 50 lines of the server log.
        """

        if kind == 'server':
            log_file = SERVER_LOG
        elif kind == 'worker':
            log_file = WORKER_LOG
        try:
            log = os.popen('tail -n %s %s' % (n, log_file)).read()
        except:
            log = "Error reading %s" % log_file

        return log
