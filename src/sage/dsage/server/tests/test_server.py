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
import tempfile
from glob import glob

from sage.dsage.database.job import Job, expand_job
from sage.dsage.database.jobdb import JobDatabaseSA as JobDatabase
from sage.dsage.database.workerdb import WorkerDatabaseSA as WorkerDatabase
from sage.dsage.database.clientdb import ClientDatabaseSA as ClientDatabase
from sage.dsage.database.db_config import init_db_sa as init_db
from sage.dsage.server.server import DSageServer

class DSageServerTestCase(unittest.TestCase):
    """
    Tests for DSageServer go here.

    """

    test_db = tempfile.NamedTemporaryFile()

    def setUp(self):
        Session = init_db(self.test_db.name)
        self.jobdb = JobDatabase(Session)
        self.workerdb = WorkerDatabase(Session)
        self.clientdb = ClientDatabase(Session)
        self.dsage_server = DSageServer(self.jobdb,
                                        self.workerdb,
                                        self.clientdb,
                                        log_level=5)
        for job in self.create_jobs(10):
            self.dsage_server.submit_job(job._reduce())

    def tearDown(self):
        self.dsage_server.jobdb._shutdown()
        self.dsage_server.jobdb.sess.close()
        from sqlalchemy.orm import clear_mappers
        clear_mappers()
        os.remove(self.test_db.name)

    def create_jobs(self, n):
        """This method creates n jobs. """

        jobs = []
        for i in range(n):
            jobs.append(Job(name='unittest', username='yqiang', code='2+2'))

        return jobs

    def testget_job(self):
        jdict = self.dsage_server.get_job()
        self.assertEquals(type(jdict), dict)
        job = expand_job(jdict)
        self.assert_(isinstance(job, Job))

    def testget_job_by_id(self):
        job = Job()
        job.code = '2+2'
        job_id = self.dsage_server.submit_job(job._reduce())
        self.assertEquals(type(job_id), str)
        self.assertEquals(len(job_id), 10)

    def testget_job_result_by_id(self):
        job = expand_job(self.dsage_server.get_job())
        job.result = 'test'
        jdict = job._reduce()
        job_id = self.dsage_server.submit_job(jdict)
        result = self.dsage_server.get_job_result_by_id(job_id)
        self.assertEquals(result, 'test')

    def testget_jobs_by_username(self):
        self.assertEquals(
                type(self.dsage_server.get_jobs_by_username('yqiang', 'new')),
                list)
        self.assertEquals(
                len(self.dsage_server.get_jobs_by_username('test', 'new')),
                0)

        job = expand_job(self.dsage_server.get_job())
        job.username = 'testing123'
        job.code = ''
        jdict = self.dsage_server.submit_job(job._reduce())
        jobs = self.dsage_server.get_jobs_by_username('testing123', 'processing')
        j = expand_job(jobs[0])
        self.assertEquals(j.username, job.username)

    def testsubmit_job(self):
        jobs = self.create_jobs(10)
        for job in jobs:
            job_id = self.dsage_server.submit_job(job._reduce())
            self.assertEquals(type(job_id), str)
            j = expand_job(self.dsage_server.get_job_by_id(job_id))
            self.assert_(isinstance(j, Job))

    def testget_all_jobs(self):
        jobs = self.dsage_server.get_all_jobs()
        self.assertEquals(len(jobs), 10)

    def testget_active_jobs(self):
        jdicts = self.dsage_server.get_all_jobs()
        for jdict in jdicts:
            job = expand_job(jdict)
            job.status = 'processing'
            jdict = job._reduce()
            job_id = self.dsage_server.submit_job(jdict)

        jdicts = self.dsage_server.get_active_jobs()
        self.assert_(len(jdicts) == 10)
        for jdict in jdicts:
            job = expand_job(jdict)
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
            self.dsage_server.submit_job(job._reduce())

        jobs = self.dsage_server.get_killed_jobs_list()
        self.assertEquals(len(jobs), 10)

        for job in jobs:
            job = expand_job(job)
            self.assert_(isinstance(job, Job))
            self.assertEquals(job.killed, True)
            self.assert_(job.update_time < datetime.datetime.now())

    def testjob_done(self):
        import time, zlib, cPickle
        job = expand_job(self.dsage_server.get_job())
        result = zlib.compress(cPickle.dumps('done'))
        output = 'done '
        completed = True
        job_id = self.dsage_server.job_done(job.job_id,
                                            output, result,
                                            completed,
                                            time.time() - time.time())
        job = expand_job(self.dsage_server.get_job_by_id(job_id))
        self.assertEquals(job.output, output)
        self.assertEquals(job.result, cPickle.loads(zlib.decompress(result)))
        self.assertEquals(job.status, 'completed')

        job = expand_job(self.dsage_server.get_job())
        result = ['testing', '123']
        output = 'testing'
        completed = False
        job_id = self.dsage_server.job_done(job.job_id, output, result,
                                            completed,
                                            time.time() - time.time())
        job = expand_job(self.dsage_server.get_job_by_id(job_id))
        self.assert_(isinstance(job.output, str))
        self.assert_(job.status != 'completed')

    def testjob_failed(self):
        job = expand_job(self.dsage_server.get_job())
        self.dsage_server.job_failed(job.job_id, 'Failure')
        job = expand_job(self.dsage_server.get_job_by_id(job.job_id))
        self.assertEquals(job.failures, 1)
        self.assertEquals(job.output, 'Failure')
        self.dsage_server.job_failed(job.job_id, 'Another Failure')
        job = expand_job(self.dsage_server.get_job_by_id(job.job_id))
        self.assertEquals(job.failures, 2)
        self.assertEquals(job.output, 'Another Failure')

    def testkill_job(self):
        job = expand_job(self.dsage_server.get_job())
        job_id = self.dsage_server.kill_job(job.job_id)
        job = expand_job(self.dsage_server.get_job_by_id(job_id))
        self.assertEquals(job.killed, True)

