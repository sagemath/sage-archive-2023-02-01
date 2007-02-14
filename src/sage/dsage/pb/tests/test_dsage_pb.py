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

import unittest
import datetime
import os
from glob import glob
from random import randint
from cPickle import dumps, loads
import zlib

from dsage.database.job import Job
from dsage.database.jobdb import JobDatabaseZODB
from dsage.pb.dsage_pb import DSage

class DSageTestCase(unittest.TestCase):
    def unpickle(self, pickled_job):
        return loads(zlib.decompress(pickled_job))

    def setUp(self):
        self.jobdb = JobDatabaseZODB(test=True)
        self.dsage = DSage(self.jobdb, log_level=5)
        for job in self.createJobs(10):
            self.dsage.jobdb.newJob(job)

    def tearDown(self):
        self.dsage.jobdb._shutdown()
        files = glob('*.db*')
        for file in files:
            os.remove(file)

    def testgetJob(self):
        pickled_job = self.dsage.getJob()
        self.assertEquals(type(pickled_job), str)
        job = self.unpickle(pickled_job)
        self.assert_(isinstance(job, Job))

    def testgetJobByID(self):
        job = self.createJobs(1)
        for j in job:
            jobID = self.jobdb.storeJob(j)
        self.assert_(isinstance(jobID, str))

    def testgetJobResultsByID(self):
        job = self.unpickle(self.dsage.getJob())
        job.result = 'test'
        job.file = ''
        id = self.dsage.submitJob(job.pickle())
        self.assert_(self.dsage.getJobResultByID(id) == 'test')

    def testgetJobsByAuthor(self):
        self.assert_(isinstance(
                     self.dsage.getJobsByAuthor('Yi Qiang',
                                                'unittest',
                                                True), list))

        self.assert_(len(self.dsage.getJobsByAuthor('test',
                                                    False,
                                                    None)) == 0)

        job = self.unpickle(self.dsage.getJob())
        job.author = 'test'
        job.file = ''
        id = self.dsage.submitJob(job.pickle())
        self.assert_(self.dsage.getJobsByAuthor('test',
                                                False,
                                                None)[0].author == 'test')

    def testsubmitJob(self):
        jobs = self.createJobs(10)
        for job in jobs:
            job.file = ""
            id = self.dsage.submitJob(job.pickle())
            self.assertEquals(type(id), str)
            j = self.unpickle(self.dsage.getJobByID(id))
            self.assert_(isinstance(j, Job))

    def testgetJobsList(self):
        jobs = self.dsage.getJobsList()
        self.assertEquals(len(jobs), 10)
        for i in xrange(len(jobs)-1):
            self.assert_(isinstance(jobs[i], Job))
            self.assert_(jobs[i].num < jobs[i+1].num)

    def testgetActiveJobs(self):
        jobs = self.dsage.getJobsList()
        for job in jobs:
            job.status = 'processing'
            id = self.dsage.submitJob(job.pickle())
        jobs = self.dsage.getActiveJobs()
        self.assert_(len(jobs) == 10)
        for job in jobs:
            self.assert_(isinstance(job, Job))
            self.assert_(job.status == 'processing')
            self.assert_(job.updated_time < datetime.datetime.now())

    def testgetActiveClientsList(self):
        pass

    def testgetKilledJobsList(self):
        jobs = self.dsage.getKilledJobsList()
        self.assertEquals(len(jobs), 0)

        jobs = self.dsage.getJobsList()
        for job in jobs:
            job.killed = True
            self.dsage.submitJob(job.pickle())

        jobs = self.dsage.getKilledJobsList()
        self.assertEquals(len(jobs), 10)

        for job in jobs:
            job = self.unpickle(job)
            self.assert_(isinstance(job, Job))
            self.assertEquals(job.killed, True)
            self.assert_(job.updated_time < datetime.datetime.now())

    def testgetNextJobID(self):
        id = self.dsage.getNextJobID()
        self.assertEquals(type(id), str)
        self.assert_(id != self.dsage.getNextJobID())

    def testjobDone(self):
        job = self.unpickle(self.dsage.getJob())
        result = 'done'
        output = 'done '
        completed = True
        id = self.dsage.jobDone(job.id, output, result, completed,
                                ('yi@test', 'no info provided'))
        job = self.unpickle(self.dsage.getJobByID(id))
        self.assertEquals(job.output, output)
        self.assertEquals(job.result, result)
        self.assertEquals(job.status, 'completed')

        job = self.unpickle(self.dsage.getJob())
        result = ['testing', '123']
        output = 'testing'
        completed = False
        id = self.dsage.jobDone(job.id, output, result, completed,
                                ('yi@test', 'no info provided'))
        job = self.unpickle(self.dsage.getJobByID(id))
        self.assert_(isinstance(job.output, str))
        self.assert_(job.status != 'completed')

    def testjobFailed(self):
        job = self.unpickle(self.dsage.getJob())
        self.dsage.jobFailed(job.id)
        job = self.unpickle(self.dsage.getJobByID(job.id))
        self.assertEquals(job.failures, 1)
        self.dsage.jobFailed(job.id)
        job = self.unpickle(self.dsage.getJobByID(job.id))
        self.assertEquals(job.failures, 2)

    def testkillJob(self):
        job = self.unpickle(self.dsage.getJob())
        reason = 'test'
        id = self.dsage.killJob(job.id, reason)
        job = self.unpickle(self.dsage.getJobByID(id))
        self.assertEquals(job.killed, True)

    def createJobs(self, n):
        """This method creates n jobs. """

        jobs = []
        for i in range(n):
            jobs.append(Job(name='unittest', author='Yi Qiang'))

        return jobs

