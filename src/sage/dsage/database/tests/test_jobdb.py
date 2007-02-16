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
import os
import datetime
from glob import glob

from dsage.database.jobdb import JobDatabaseZODB, DatabasePruner
from dsage.database.job import Job

class JobDatabaseZODBTestCase(unittest.TestCase):
    def setUp(self):
        self.jobdb = JobDatabaseZODB(test=True)

        jobs = self.createJobs(10)
        for job in jobs:
            self.jobdb.newJob(job)

    def tearDown(self):
        self.jobdb._shutdown()
        files = glob('*.db*')
        for file in files:
            os.remove(file)

    def testhasJob(self):
        job = self.createJobs(1)
        jobID = self.jobdb.storeJob(job[0])
        self.assertEquals(self.jobdb.hasJob(jobID), True)
        self.assert_(self.jobdb.hasJob('GLOBALS') > 0)

    def testgetJobID(self):
        """Attempt to increment Job ID by one. """

        jobID = self.jobdb.getJobID()
        jobID1 = self.jobdb.getJobID()

        self.assert_(int(jobID[11:]) < int(jobID1[11:]))

    def testgetJob(self):
        job = self.jobdb.getJob()
        self.assert_(isinstance(job, Job) and job.status != 'completed')

    def testremoveJob(self):
        self.assertRaises(KeyError, self.jobdb.removeJob, 'not_a_job_id')

        job = self.jobdb.getJob()
        self.assertEquals(type(job), Job)
        job_id = job.id
        self.assertEquals(self.jobdb.removeJob(job_id), job_id)
        self.assertRaises(KeyError, self.jobdb.removeJob, job_id)

    def teststoreJob(self):
        jobs = self.createJobs(10)
        for job in jobs:
            job_id = self.jobdb.newJob(job)
            self.assertEquals(job, self.jobdb.getJobByID(job_id))
            self.assert_(self.jobdb.getJobByID(job_id).updated_time <
                         datetime.datetime.now())

    def testgetJobsByAuthor(self):
        jobs = self.jobdb.getJobsByAuthor('Yi Qiang', True)
        self.assert_(len(jobs) == 0)

        jobs = self.jobdb.getJobsByAuthor('Yi Qiang', False, 'unittest')
        self.assert_(len(jobs) > 0)

    def testgetActiveJobs(self):
        jobs = self.jobdb.getActiveJobs()
        self.assert_(len(jobs) == 0)

        jobs = self.jobdb.getJobsList()
        for job in jobs:
            job.status = 'processing'
            self.jobdb.storeJob(job)

        jobs = self.jobdb.getActiveJobs()
        self.assert_(len(jobs) == 10)

    def testgetJobsList(self):
        jobs = self.jobdb.getJobsList()

        for i in xrange(len(jobs)-1):
            self.assert_(jobs[i].num < jobs[i+1].num)

    def createJobs(self, n):
        """This method creates n jobs. """

        jobs = []
        for i in range(n):
            jobs.append(Job(name='unittest', author='Yi Qiang'))

        return jobs

class DatabasePrunerTestCase(unittest.TestCase):
    def setUp(self):
        self.jobdb = JobDatabaseZODB(test=True)
        self.pruner = DatabasePruner(self.jobdb)
        jobs = self.createJobs(10)
        for job in jobs:
            self.jobdb.newJob(job)


    def tearDown(self):
        self.jobdb._shutdown()
        files = glob('*.db*')
        for file in files:
            os.remove(file)

    def testcleanOldJobs(self):
        jobs = self.jobdb.getJobsList()
        for job in jobs:
            job.updated_time -= datetime.timedelta(10)
            # directly accessing the database because storeJob
            # automatically updates the updated_time
            self.jobdb.jobdb[job.id] = job

        self.pruner.cleanOldJobs()

        jobs = self.jobdb.getJobsList()
        self.assertEquals(len(jobs), 0)

    def testcleanFailedJobs(self):
        jobs = self.jobdb.getJobsList()
        for job in jobs:
            job.failures += 20
            self.jobdb.storeJob(job)

        self.pruner.cleanFailedJobs()
        jobs = self.jobdb.getJobsList()
        self.assertEquals(len(jobs), 0)

    def createJobs(self, n):
        """This method creates n jobs. """

        jobs = []
        for i in range(n):
            jobs.append(Job(name='unittest', author='Yi Qiang'))

        return jobs

if __name__ == '__main__':
    unittest.main()
