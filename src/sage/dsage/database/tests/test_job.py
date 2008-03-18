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

from sage.dsage.database.job import Job
from sage.dsage.misc.misc import random_string

class JobTestCase(unittest.TestCase):
    def testcreate_job(self):
        job = Job()
        self.assert_(isinstance(job, Job))
        job_id = random_string(10)
        job = Job(job_id=job_id, name='test', code='test', username='test',
                  kind='test')
        self.assert_(isinstance(job, Job))

    def testjob_id(self):
        job = Job()
        self.assertRaises(TypeError, job.job_id, 5)

    def testjobFile(self):
        job = Job()
        self.assertRaises(TypeError, job.code, 1)

    def testjobCreationTime(self):
        job = Job()
        self.assert_(isinstance(job.creation_time, datetime.datetime))
        self.assertRaises(TypeError, job.creation_time, 'test')

    def testjobFinishTime(self):
        job = Job()
        self.assert_(job.finish_time == None)
        self.assertRaises(TypeError, job.finish_time, 'test')

    def testjobUpdateTime(self):
        job = Job()
        self.assertEquals(job.update_time, None)
        self.assertRaises(TypeError, job.update_time, 'test')

    def testjobStatus(self):
        job = Job()
        self.assertRaises(TypeError, job.status, 'test')

    def testjobFailures(self):
        job = Job()
        self.assertRaises(TypeError, job.failures, 'test')

    def testjobResult(self):
        job = Job()
        job.result = 'woot'
        self.assertEquals(job.result, 'woot')

    def testjobKilled(self):
        job = Job()
        self.assertRaises(TypeError, job.killed, 'test')


