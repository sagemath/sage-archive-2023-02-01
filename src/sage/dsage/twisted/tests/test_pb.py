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

from sage.dsage.database.job import Job
from sage.dsage.database.jobdb import JobDatabaseZODB
from sage.dsage.server.server import DSageServer

class DSageTestCase(unittest.TestCase):
    def unpickle(self, pickled_job):
        return loads(zlib.decompress(pickled_job))

    def setUp(self):
        self.jobdb = JobDatabaseZODB(test=True)
        self.dsage_server = DSageServer(self.jobdb, log_level=5)
        for job in self.create_jobs(10):
            self.dsage_server.jobdb.new_job(job)

    def tearDown(self):
        self.dsage_server.jobdb._shutdown()
        files = glob('*.db*')
        for file in files:
            os.remove(file)

    def testget_job(self):
        pickled_job = self.dsage_server.get_job()
        self.assertEquals(type(pickled_job), str)
        job = self.unpickle(pickled_job)
        self.assert_(isinstance(job, Job))

    def testget_job_by_id(self):
        job = self.create_jobs(1)
        for j in job:
            jobID = self.jobdb.store_job(j)
        self.assert_(isinstance(jobID, str))

    def testget_job_resultsByID(self):
        job = self.unpickle(self.dsage_server.get_job())
        job.result = 'test'
        job.code = ''
        id = self.dsage_server.submit_job(job.pickle())
        self.assert_(self.dsage_server.get_job_result_by_id(id) == 'test')

    def testget_jobs_by_user_id(self):
        self.assert_(isinstance(
                     self.dsage_server.get_jobs_by_user_id('Yi Qiang',
                                                'unittest',
                                                True), list))

        self.assert_(len(self.dsage_server.get_jobs_by_user_id('test',
                                                    False,
                                                    None)) == 0)

        job = self.unpickle(self.dsage_server.get_job())
        job.user_id = 'test'
        job.code = ''
        id = self.dsage_server.submit_job(job.pickle())
        self.assert_(self.dsage_server.get_jobs_by_user_id('test',
                                                False,
                                                None)[0].user_id == 'test')

    def testsubmit_job(self):
        jobs = self.create_jobs(10)
        for job in jobs:
            job.code = ""
            id = self.dsage_server.submit_job(job.pickle())
            self.assertEquals(type(id), str)
            j = self.unpickle(self.dsage_server.get_job_by_id(id))
            self.assert_(isinstance(j, Job))

    def testget_jobs_list(self):
        jobs = self.dsage_server.get_jobs_list()
        self.assertEquals(len(jobs), 10)
        for i in xrange(len(jobs)-1):
            self.assert_(isinstance(jobs[i], Job))
            self.assert_(jobs[i].num < jobs[i+1].num)

    def testget_active_jobs(self):
        jobs = self.dsage_server.get_jobs_list()
        for job in jobs:
            job.status = 'processing'
            id = self.dsage_server.submit_job(job.pickle())
        jobs = self.dsage_server.get_active_jobs()
        self.assert_(len(jobs) == 10)
        for job in jobs:
            self.assert_(isinstance(job, Job))
            self.assert_(job.status == 'processing')
            self.assert_(job.update_time < datetime.datetime.now())

    def testget_active_clients_list(self):
        pass

    def testget_killed_jobs_list(self):
        jobs = self.dsage_server.get_killed_jobs_list()
        self.assertEquals(len(jobs), 0)

        jobs = self.dsage_server.get_jobs_list()
        for job in jobs:
            job.killed = True
            self.dsage_server.submit_job(job.pickle())

        jobs = self.dsage_server.get_killed_jobs_list()
        self.assertEquals(len(jobs), 10)

        for job in jobs:
            job = self.unpickle(job)
            self.assert_(isinstance(job, Job))
            self.assertEquals(job.killed, True)
            self.assert_(job.update_time < datetime.datetime.now())

    def testget_next_job_id(self):
        id = self.dsage_server.get_next_job_id()
        self.assertEquals(type(id), str)
        self.assert_(id != self.dsage_server.get_next_job_id())

    def testjob_done(self):
        job = self.unpickle(self.dsage_server.get_job())
        result = 'done'
        output = 'done '
        completed = True
        id = self.dsage_server.job_done(job.id, output, result, completed,
                                ('yi@test', 'no info provided'))
        job = self.unpickle(self.dsage_server.get_job_by_id(id))
        self.assertEquals(job.output, output)
        self.assertEquals(job.result, result)
        self.assertEquals(job.status, 'completed')

        job = self.unpickle(self.dsage_server.get_job())
        result = ['testing', '123']
        output = 'testing'
        completed = False
        id = self.dsage_server.job_done(job.id, output, result, completed,
                                ('yi@test', 'no info provided'))
        job = self.unpickle(self.dsage_server.get_job_by_id(id))
        self.assert_(isinstance(job.output, str))
        self.assert_(job.status != 'completed')

    def testjob_failed(self):
        job = self.unpickle(self.dsage_server.get_job())
        self.dsage_server.job_failed(job.id)
        job = self.unpickle(self.dsage_server.get_job_by_id(job.id))
        self.assertEquals(job.failures, 1)
        self.dsage_server.job_failed(job.id)
        job = self.unpickle(self.dsage_server.get_job_by_id(job.id))
        self.assertEquals(job.failures, 2)

    def testkill_job(self):
        job = self.unpickle(self.dsage_server.get_job())
        reason = 'test'
        id = self.dsage_server.kill_job(job.id, reason)
        job = self.unpickle(self.dsage_server.get_job_by_id(id))
        self.assertEquals(job.killed, True)

    def create_jobs(self, n):
        """This method creates n jobs. """

        jobs = []
        for i in range(n):
            jobs.append(Job(name='unittest', user_id='Yi Qiang'))

        return jobs

