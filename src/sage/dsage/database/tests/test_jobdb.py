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

import unittest
import os
import datetime
from glob import glob

from sage.dsage.database.jobdb import JobDatabaseZODB, JobDatabaseSQLite
from sage.dsage.database.jobdb import DatabasePruner
from sage.dsage.database.job import Job

class JobDatabaseZODBTestCase(unittest.TestCase):
    def setUp(self):
        self.jobdb = JobDatabaseZODB(test=True)

        jobs = self.create_jobs(10)
        for job in jobs:
            self.jobdb.new_job(job)

    def tearDown(self):
        self.jobdb._shutdown()
        files = glob('*.db*')
        for file in files:
            os.remove(file)

    def testhas_job(self):
        job = self.create_jobs(1)
        job_id = self.jobdb.store_job(job[0])
        self.assertEquals(self.jobdb.has_job(job_id), True)
        self.assert_(self.jobdb.has_job('GLOBALS') > 0)

    def testget_next_job_id(self):
        """Attempt to increment Job ID by one. """

        job_id = self.jobdb.get_next_job_id()
        job_id1 = self.jobdb.get_next_job_id()

        self.assert_(int(job_id[11:]) < int(job_id1[11:]))

    def testget_job(self):
        job = self.jobdb.get_job()
        self.assert_(isinstance(job, Job) and job.status != 'completed')

    def testremove_job(self):
        self.assertRaises(KeyError, self.jobdb.remove_job, 'not_a_job_id')

        job = self.jobdb.get_job()
        self.assertEquals(type(job), Job)
        job_id = job.job_id
        self.assertEquals(self.jobdb.remove_job(job_id), job_id)
        self.assertRaises(KeyError, self.jobdb.remove_job, job_id)

    def teststore_job(self):
        jobs = self.create_jobs(10)
        for job in jobs:
            job_id = self.jobdb.new_job(job)
            self.assertEquals(job, self.jobdb.get_job_by_id(job_id))
            self.assert_(self.jobdb.get_job_by_id(job_id).update_time <=
                         datetime.datetime.now())

    def testget_jobs_by_username(self):
        jobs = self.jobdb.get_jobs_by_username('Yi Qiang', True)
        self.assert_(len(jobs) == 0)

        jobs = self.jobdb.get_jobs_by_username('Yi Qiang', False, 'unittest')
        self.assert_(len(jobs) > 0)

    def testget_active_jobs(self):
        jobs = self.jobdb.get_active_jobs()
        self.assert_(len(jobs) == 0)

        jobs = self.jobdb.get_jobs_list()
        for job in jobs:
            job.status = 'processing'
            self.jobdb.store_job(job)

        jobs = self.jobdb.get_active_jobs()
        self.assert_(len(jobs) == 10)

    def testget_jobs_list(self):
        jobs = self.jobdb.get_jobs_list()
        self.assertEquals(len(jobs), 10)

    def create_jobs(self, n):
        """This method creates n jobs. """

        jobs = []
        for i in range(n):
            jobs.append(Job(name='unittest', username='Yi Qiang'))

        return jobs

class JobDatabaseSQLiteTestCase(unittest.TestCase):
    """
    Unit tests for the SQLite based JobDatabase go here.

    """

    def setUp(self):
        self.jobdb = JobDatabaseSQLite(test=True)

    def tearDown(self):
        query = """DELETE FROM jobs"""
        cur = self.jobdb.con.cursor()
        cur.execute(query)
        self.jobdb._shutdown()

    def testget_job(self):
        job = Job()
        job.status = 'new'
        job.killed = False
        job_id = self.jobdb.store_job(job.reduce())
        self.assertEquals(job_id, self.jobdb.get_job()['job_id'])

    def teststore_job(self):
        job = Job()
        self.assert_(isinstance(self.jobdb.store_job(job.reduce()), str))

    def testget_job_by_id(self):
        job = Job()
        job_id = self.jobdb.store_job(job.reduce())
        self.assert_(self.jobdb.get_job_by_id(job_id) is not None)

    def testhas_job(self):
        job = Job()
        job_id = self.jobdb.store_job(job.reduce())
        self.assertEquals(self.jobdb.has_job(job_id), True)

    def testcreate_jdict(self):
        job = Job()
        job_id = self.jobdb.store_job(job.reduce())
        jdict = self.jobdb.get_job_by_id(job_id)
        self.assert_(isinstance(jdict, dict))

    def testget_killed_jobs_list(self):
        job = Job()
        job_id = self.jobdb.store_job(job.reduce())
        jdict = self.jobdb.get_job_by_id(job_id)
        self.jobdb.set_killed(jdict['job_id'], killed=True)
        self.assertEquals(len(self.jobdb.get_killed_jobs_list()), 1)

class DatabasePrunerTestCase(unittest.TestCase):
    def setUp(self):
        self.jobdb = JobDatabaseZODB(test=True)
        self.pruner = DatabasePruner(self.jobdb)
        jobs = self.create_jobs(10)
        for job in jobs:
            self.jobdb.new_job(job)


    def tearDown(self):
        self.jobdb._shutdown()
        files = glob('*.db*')
        for file in files:
            os.remove(file)

    def testclean_old_jobs(self):
        jobs = self.jobdb.get_jobs_list()
        for job in jobs:
            job.update_time -= datetime.timedelta(10)
            # directly accessing the database because store_job
            # automatically updates the update_time
            self.jobdb.jobdb[job.job_id] = job

        self.pruner.clean_old_jobs()

        jobs = self.jobdb.get_jobs_list()
        self.assertEquals(len(jobs), 0)

    def testclean_failed_jobs(self):
        jobs = self.jobdb.get_jobs_list()
        for job in jobs:
            job.failures += 20
            self.jobdb.store_job(job)

        self.pruner.clean_failed_jobs()
        jobs = self.jobdb.get_jobs_list()
        self.assertEquals(len(jobs), 0)

    def create_jobs(self, n):
        """This method creates n jobs. """

        jobs = []
        for i in range(n):
            jobs.append(Job(name='unittest', username='Yi Qiang'))

        return jobs

if __name__ == '__main__':
    unittest.main()
