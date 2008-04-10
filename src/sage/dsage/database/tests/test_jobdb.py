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
import tempfile

from sage.dsage.database.jobdb import JobDatabaseSQLite, JobDatabaseSA
from sage.dsage.database.job import Job
from sage.dsage.database.db_config import create_schema, init_db, init_db_sa

class JobDatabaseSATestCase(unittest.TestCase):
    test_db = tempfile.NamedTemporaryFile()

    def setUp(self):
        Session = init_db_sa(self.test_db.name)
        self.jobdb = JobDatabaseSA(Session)

    def tearDown(self):
        from sqlalchemy.orm import clear_mappers
        self.jobdb.sess.close()
        clear_mappers()
        os.remove(self.test_db.name)

    def testget_job(self):
        job = Job()
        job.status = 'new'
        job.killed = False
        job_id = self.jobdb.store_jdict(job._reduce())
        self.assertEquals(job_id, self.jobdb.get_job().job_id)

    def teststore_jdict(self):
        job = Job()
        self.assert_(isinstance(self.jobdb.store_jdict(job._reduce()), str))


class JobDatabaseSQLiteTestCase(unittest.TestCase):
    """
    Unit tests for the SQLite based JobDatabase go here.

    """

    test_db = tempfile.NamedTemporaryFile()

    def setUp(self):
        db_conn = init_db(self.test_db.name)
        create_schema(db_conn)
        self.jobdb = JobDatabaseSQLite(db_conn)

    def tearDown(self):
        query = """DELETE FROM jobs"""
        cur = self.jobdb.con.cursor()
        cur.execute(query)
        self.jobdb._shutdown()
        os.remove(self.test_db.name)

    def testget_job(self):
        job = Job()
        job.status = 'new'
        job.killed = False
        job_id = self.jobdb.store_jdict(job._reduce())
        self.assertEquals(job_id, self.jobdb.get_job()['job_id'])

    def teststore_jdict(self):
        job = Job()
        self.assert_(isinstance(self.jobdb.store_jdict(job._reduce()), str))

    def testget_job_by_id(self):
        job = Job()
        job_id = self.jobdb.store_jdict(job._reduce())
        self.assert_(self.jobdb.get_job_by_id(job_id) is not None)

    def testhas_job(self):
        job = Job()
        job_id = self.jobdb.store_jdict(job._reduce())
        self.assertEquals(self.jobdb.has_job(job_id), True)

    def testcreate_jdict(self):
        job = Job()
        job_id = self.jobdb.store_jdict(job._reduce())
        jdict = self.jobdb.get_job_by_id(job_id)
        self.assert_(isinstance(jdict, dict))

    def testget_killed_jobs_list(self):
        job = Job()
        job_id = self.jobdb.store_jdict(job._reduce())
        jdict = self.jobdb.get_job_by_id(job_id)
        self.jobdb.set_killed(jdict['job_id'], killed=True)
        self.assertEquals(len(self.jobdb.get_killed_jobs_list()), 1)

if __name__ == '__main__':
    unittest.main()
