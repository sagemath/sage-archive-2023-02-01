############################################################################
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
############################################################################

import datetime
import copy
import cPickle
import zlib

from twisted.internet import reactor, task

from sage.dsage.database.job import Job
from sage.dsage.interface.dsage_interface import JobWrapper, blockingJobWrapper
from sage.dsage.twisted.misc import blockingCallFromThread

class DistributedFunction(object):
    r"""
    Parent class for all classes that wish to use Distributed SAGE.

    Parameters:
        dsage -- a DSage() object

    """

    def __init__(self, dsage):
        self.DSage = dsage
        self.result = None
        self.outstanding_jobs = [] # unsubmitted jobs
        self.waiting_jobs = [] # unprocessed jobs
        self.processed_jobs = [] # processed jobs
        self.start_time = None
        self.end_time= None
        self.name = None
        self.job_files = []
        # self.checker_task = None
        self._done = False

    def __getstate__(self):
        d = copy(self.__dict__)
        d['DSage'] = None
        d['checker_task'] = None
        return d

    def get_done(self):
        return self._done
    def set_done(self, value):
        if value:
            self.end_time = datetime.datetime.now()
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
        r"""
        Saves your distributed job to a sobj.

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

    def restore(self, dsage):
        if dsage.remoteobj is None:
            # XXX This is a hack because dsage.remoteobj is not set yet
            from twisted.internet import reactor
            reactor.callLater(1.0, self.restore, dsage)
            return
        self.DSage = dsage
        for job in self.waiting_jobs:
            job.remoteobj = self.DSage.remoteobj
        if not len(self.outstanding_jobs) == 0:
            self.submit_jobs(self.name)
        self.checker_task = task.LoopingCall(self.check_results)
        self.checker_task.start(2.0, now=True)

    def submit_job(self, job, job_name='job', async=False):
        if async:
            if isinstance(job, Job):
                self.waiting_jobs.append(self.DSage.send_job(job,
                                                             async=True))
            else:
                self.waiting_jobs.append(self.DSage.eval(job,
                                                         job_name=job_name,
                                                         async=True))
        else:
            if isinstance(job, Job):
                self.waiting_jobs.append(self.DSage.send_job(job))
            else:
                self.waiting_jobs.append(self.DSage.eval(job,
                                                         job_name=job_name))

    def submit_jobs(self, job_name='job', async=False):
        for job in self.outstanding_jobs:
           self.submit_job(job, job_name, async)
        self.outstanding_jobs = []

    def start(self):
        self.start_time = datetime.datetime.now()
        reactor.callFromThread(self.submit_jobs, self.name, async=True)
        self.checker_task = blockingCallFromThread(task.LoopingCall,
                                                   self.check_results)
        reactor.callFromThread(self.checker_task.start,
                               1.0, now=True)
    def check_results(self):
        for wrapped_job in self.waiting_jobs:
            if isinstance(wrapped_job, JobWrapper):
                wrapped_job.getJob()
            else:
                wrapped_job.async_getJob()
            if wrapped_job.status == 'completed':
                self.waiting_jobs.remove(wrapped_job)
                self.process_result(wrapped_job)
                self.processed_jobs.append(wrapped_job)
        if self.done:
            # kill the jobs in the waiting queue
            for wrapped_job in self.waiting_jobs:
                if isinstance(wrapped_job, JobWrapper):
                    wrapped_job.kill()
                else:
                    wrapped_job.async_kill()
            self.waiting_jobs = []
            reactor.callFromThread(self.checker_task.stop)

        self.done = len(self.waiting_jobs) == 0

class DistributedFunctionTest(DistributedFunction):
    def __init__(self, DSage, n, name='DistributedFunctionTest'):
        DistributedFunction.__init__(self, DSage)
        self.n = n
        self.name = name
        self.result = 0
        self.results = []
        self.outstanding_jobs = ["print %s"%i for i in range(1,n+1)]

    def process_result(self, job):
        self.result += int(job.output)