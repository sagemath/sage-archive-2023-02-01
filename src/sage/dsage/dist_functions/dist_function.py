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
#
##############################################################################
import time
import datetime
import cPickle
import zlib

from sage.dsage.database.job import Job
from sage.dsage.interface.dsage_interface import (JobWrapper,
                                                  BlockingJobWrapper,
                                                  blockingCallFromThread)

class DistributedFunction(object):
    """
    Parent class for all classes that wish to use Distributed SAGE.

    Parameters:
        dsage -- a DSage() object

    """

    def __init__(self, dsage):
        self._dsage = dsage
        self.result = None
        self.outstanding_jobs = [] # unsubmitted jobs
        self.waiting_jobs = [] # unprocessed jobs
        self.processed_jobs = [] # processed jobs
        self.start_time = None
        self.end_time= None
        self.name = None
        self.job_files = []
        self._done = False

    def __getstate__(self):
        import copy
        d = copy.copy(self.__dict__)
        d['DSage'] = None
        d['checker_task'] = None

        return d

    def get_done(self):
        return self._done
    def set_done(self, value):
        if value:
            self.end_time = datetime.datetime.now()
            # self.stop()
        self._done = value
    done = property(fget=get_done, fset=set_done, fdel=None,
                    doc="Set the jobs status")


    def get_time(self):
        if self.end_time == None:
            return "End time not available"
        return str(self.end_time - self.start_time)
    time = property(fget=get_time, fset=None, fdel=None,
                     doc="Time it took for job to complete")


    def save(self, filename=None, compress=True):
        """
        Saves your distributed job to disk.

        """

        if filename is None:
            filename = str(self.name)
        if not filename.endswith('.sobj'):
            filename += '.sobj'
        s = cPickle.dumps(self, 2)
        if compress:
            s = zlib.compress(s)
        f = open(filename, 'wb')
        f.write(s)
        f.close()

        return filename


    def stop(self, verbose=True):
        """
        Ends the current DistributedFunction, kills all waiting jobs.

        """

        for job in self.waiting_jobs:
            job.kill()
        self.done = True
        if verbose:
            print 'All waiting jobs have been killed. This job is no more.'


    def restore(self, dsage):
        """
        Reloads a distributed job from disk.

        """
        from twisted.internet import reactor
        from twisted.internet import task
        if dsage.remoteobj is None:
            # XXX This is a hack because dsage.remoteobj is not set yet
            self.reactor.callLater(0.5, self.restore, dsage)
            return
        self._dsage = dsage
        for job in self.waiting_jobs:
            job.remoteobj = self._dsage.remoteobj
        if not len(self.outstanding_jobs) == 0:
            self.submit_jobs(self.name)

        self.checker_task = task.LoopingCall(self.check_waiting_jobs)
        self.reactor.callFromThread(self.checker_task.start, 5.0, now=True)


    def submit_job(self, job, job_name='job', async=True):
        """
        Submits a job to the server.

        """

        job.username = self._dsage.username
        if async:
            if isinstance(job, Job):
                self.waiting_jobs.append(self._dsage.send_job(job,
                                         async=True))
            else:
                self.waiting_jobs.append(self._dsage.eval(job,
                                                         job_name=job_name,
                                                         async=True))
        else:
            if isinstance(job, Job):
                self.waiting_jobs.append(self._dsage.send_job(job,
                                                             async=False))
            else:
                self.waiting_jobs.append(self._dsage.eval(job,
                                                         job_name=job_name))


    def submit_jobs(self, job_name='job', async=True):
        """
        Repeatedly calls submit_job until we have no more jobs in
        outstanding_jobs

        """

        for job in self.outstanding_jobs:
            if isinstance(job, Job):
               self.submit_job(job, job_name, async)
        self.outstanding_jobs = []

    def wait(self, timeout=None):
        """
        Blocks until the job is completed.

        """

        import signal
        if timeout == None:
            while not self.done:
                    time.sleep(0.5)
        else:
            def handler(signum, frame):
                raise RuntimeError('Maximum wait time exceeded.')
            signal.signal(signal.SIGALRM, handler)
            signal.alarm(timeout)
            while not self.done:
                time.sleep(0.5)
            signal.alarm(0)

    def start(self):
        """
        Starts the Distributed Function. It will submit all jobs in the
        outstanding_jobs queue and also start a checker tasks that polls for
        new information about jobs in the waiting_jobs queue

        """

        from twisted.internet import reactor, task
        self.reactor = reactor
        if self._dsage is None:
            print 'Error: Not connected to a DSage server.'
            return
        self.start_time = datetime.datetime.now()
        self.reactor.callFromThread(self.submit_jobs, self.name, async=True)
        self.checker_task = blockingCallFromThread(self.reactor,
                                                   task.LoopingCall,
                                                   self.check_waiting_jobs)
        self.reactor.callFromThread(self.checker_task.start, 5.0, now=True)


    def process_result(self):
        """
        Any class subclassing DistributedFunction should implement this
        method.

        """

        pass


    def check_waiting_jobs(self):
        """
        Checks the status of jobs in the waiting queue.

        """

        from twisted.internet import reactor
        from twisted.spread import pb
        for wrapped_job in self.waiting_jobs:
            if wrapped_job.killed == True:
                self.waiting_jobs.remove(wrapped_job)
                continue
            if isinstance(wrapped_job, JobWrapper):
                try:
                    wrapped_job.get_job()
                except pb.DeadReferenceError, msg:
                    if self.checker_task.running:
                        self.reactor.callFromThread(self.checker_task.stop)
                    break
            elif isinstance(wrapped_job, BlockingJobWrapper):
                wrapped_job.async_get_job()
            if wrapped_job.status == 'completed':
                self.waiting_jobs.remove(wrapped_job)
                self.process_result(wrapped_job)
                self.processed_jobs.append(wrapped_job)
        if self.done:
            for wrapped_job in self.waiting_jobs:
                if isinstance(wrapped_job, JobWrapper):
                    wrapped_job.kill()
                elif isinstance(wrapped_job, BlockingJobWrapper):
                    wrapped_job.async_kill()
            self.waiting_jobs = []
            if self.checker_task.running:
                self.reactor.callFromThread(self.checker_task.stop)


class DistributedFunctionTest(DistributedFunction):
    """
    This is a very simple DistributedFunction.
    Only for educational purposes.

    """

    def __init__(self, dsage, n, name='DistributedFunctionTest'):
        DistributedFunction.__init__(self, dsage)
        self.n = n
        self.name = name
        self.result = 0
        self.results = []
        self.code = """DSAGE_RESULT=%s"""
        self.outstanding_jobs = [Job(code=self.code % i, username='yqiang')
                                 for i in range(1, n+1)]
        self.start()

    def process_result(self, job):
        self.done = len(self.waiting_jobs) == 0
        self.result += (job.result)
