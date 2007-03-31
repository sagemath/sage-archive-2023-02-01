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

from sage.dsage.database.job import Job, expand_job
from sage.dsage.database.jobdb import JobDatabaseSQLite
from sage.dsage.database.monitordb import MonitorDatabase
from sage.dsage.database.clientdb import ClientDatabase
from sage.dsage.server.server import DSageServer

class DSageServerTestCase(unittest.TestCase):
    r"""
    Tests for DSageServer go here.

    """

    def setUp(self):
        self.jobdb = JobDatabaseSQLite(test=True)
        self.monitordb = MonitorDatabase(test=True)
        self.clientdb = ClientDatabase(test=True)
        self.dsage_server = DSageServer(self.jobdb,
                                        self.monitordb,
                                        self.clientdb,
                                        log_level=5)
        for job in self.create_jobs(10):
            self.dsage_server.submit_job(job.reduce())

    def tearDown(self):
        self.dsage_server.jobdb._shutdown()
        files = glob('*.db*')
        for file in files:
            os.remove(file)

    def testget_job(self):
        jdict = self.dsage_server.get_job()
        self.assertEquals(type(jdict), dict)
        job = expand_job(jdict)
        self.assert_(isinstance(job, Job))

    def testget_job_by_id(self):
        job = Job()
        job.code = '2+2'
        jdict = self.dsage_server.submit_job(job.reduce())
        self.assertEquals(type(jdict['job_id']), str)

    def testget_job_result_by_id(self):
        job = self.create_jobs(1)[0]
        job = expand_job(self.dsage_server.get_job())
        job.result = 'test'
        jdict = self.dsage_server.submit_job(job.reduce())
        self.assertEquals(
                    self.dsage_server.get_job_result_by_id(jdict['job_id']),
                    'test')

    def testget_jobs_by_username(self):
        self.assertEquals(
                type(self.dsage_server.get_jobs_by_username('yqiang')),
                list)
        self.assertEquals(
                len(self.dsage_server.get_jobs_by_username('test')),
                0)

        job = expand_job(self.dsage_server.get_job())
        job.username = 'testing123'
        job.code = ''
        jdict = self.dsage_server.submit_job(job.reduce())
        j = expand_job(self.dsage_server.get_jobs_by_username('testing123')[0])
        self.assertEquals(j.username, job.username)

    def testsubmit_job(self):
        jobs = self.create_jobs(10)
        for job in jobs:
            jdict = self.dsage_server.submit_job(job.reduce())
            self.assertEquals(type(jdict), dict)
            j = expand_job(self.dsage_server.get_job_by_id(jdict['job_id']))
            self.assert_(isinstance(j, Job))

    def testget_all_jobs(self):
        jobs = self.dsage_server.get_all_jobs()
        self.assertEquals(len(jobs), 10)

    def testget_active_jobs(self):
        jobs = self.dsage_server.get_all_jobs()
        for job in jobs:
            job = expand_job(job)
            job.status = 'processing'
            jdict = self.dsage_server.submit_job(job.reduce())
        jobs = self.dsage_server.get_active_jobs()
        self.assert_(len(jobs) == 10)
        for job in jobs:
            job = expand_job(job)
            self.assert_(isinstance(job, Job))
            self.assert_(job.status == 'processing')
            self.assert_(job.update_time < datetime.datetime.now())

    def testget_active_clients_list(self):
        clients = self.dsage_server.get_client_list()
        self.assertEquals(len(clients), 0)
        self.assertEquals(type(clients), list)

    def testget_killed_jobs_list(self):
        jobs = self.dsage_server.get_killed_jobs_list()
        self.assertEquals(len(jobs), 0)

        jobs = self.dsage_server.get_all_jobs()
        for job in jobs:
            job = expand_job(job)
            job.killed = True
            self.dsage_server.submit_job(job.reduce())

        jobs = self.dsage_server.get_killed_jobs_list()
        self.assertEquals(len(jobs), 10)

        for job in jobs:
            job = expand_job(job)
            self.assert_(isinstance(job, Job))
            self.assertEquals(job.killed, True)
            self.assert_(job.update_time < datetime.datetime.now())

    def testjob_done(self):
        job = expand_job(self.dsage_server.get_job())
        result = 'done'
        output = 'done '
        completed = True
        jdict = self.dsage_server.job_done(job.id, output, result, completed,
                                ('yi@test', 'no info provided'))
        job = expand_job(self.dsage_server.get_job_by_id(jdict['job_id']))
        self.assertEquals(job.output, output)
        self.assertEquals(job.result, result)
        self.assertEquals(job.status, 'completed')

        job = expand_job(self.dsage_server.get_job())
        result = ['testing', '123']
        output = 'testing'
        completed = False
        jdict = self.dsage_server.job_done(job.id, output, result, completed,
                                ('yi@test', 'no info provided'))
        job = expand_job(self.dsage_server.get_job_by_id(jdict['job_id']))
        self.assert_(isinstance(job.output, str))
        self.assert_(job.status != 'completed')

    def testjob_failed(self):
        job = expand_job(self.dsage_server.get_job())
        self.dsage_server.job_failed(job.id)
        job = expand_job(self.dsage_server.get_job_by_id(job.id))
        self.assertEquals(job.failures, 1)
        self.dsage_server.job_failed(job.id)
        job = expand_job(self.dsage_server.get_job_by_id(job.id))
        self.assertEquals(job.failures, 2)

    def testkill_job(self):
        job = expand_job(self.dsage_server.get_job())
        reason = 'test'
        id = self.dsage_server.kill_job(job.id, reason)
        job = expand_job(self.dsage_server.get_job_by_id(id))
        self.assertEquals(job.killed, True)

    def create_jobs(self, n):
        """This method creates n jobs. """

        jobs = []
        for i in range(n):
            jobs.append(Job(name='unittest', username='yqiang', code='2+2'))

        return jobs

