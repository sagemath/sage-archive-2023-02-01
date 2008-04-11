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
import os
import random
import tempfile
from glob import glob

from twisted.trial import unittest
from twisted.spread import pb
from twisted.internet import reactor, defer
from twisted.cred import portal, credentials
from twisted.conch.ssh import keys
from twisted.python import log

from sage.dsage.twisted.pb import Realm
from sage.dsage.server.server import DSageServer
from sage.dsage.twisted.pb import _SSHKeyPortalRoot
from sage.dsage.twisted.pb import ClientFactory
from sage.dsage.twisted.pubkeyauth import PublicKeyCredentialsCheckerDB
from sage.dsage.database.jobdb import JobDatabaseSA as JobDatabase
from sage.dsage.database.workerdb import WorkerDatabaseSA as WorkerDatabase
from sage.dsage.database.clientdb import ClientDatabaseSA as ClientDatabase
from sage.dsage.database.db_config import init_db_sa as init_db
from sage.dsage.database.job import Job
from sage.dsage.errors.exceptions import BadJobError
from sage.dsage.misc.hostinfo import HostInfo
from sage.dsage.misc.misc import gen_uuid
from sage.dsage.twisted.tests.test_pubkeyauth import TEST_PUB_KEY
from sage.dsage.twisted.tests.test_pubkeyauth import TEST_PRIV_KEY

hf = HostInfo().host_info
hf['uuid'] = gen_uuid()
hf['workers'] = 2

RANDOM_DATA =  ''.join([chr(i) for i in [random.randint(65, 123) for n in
                       range(500)]])

class RemoteTests(unittest.TestCase):
    test_db = tempfile.NamedTemporaryFile()

    def _login(self, *args, **kwargs):
        d = self.factory.login(self.creds, None)

    def setUp(self):
        from sqlalchemy.orm import clear_mappers
        clear_mappers()

        Session = init_db(self.test_db.name)
        self.jobdb = JobDatabase(Session)
        self.workerdb = WorkerDatabase(Session)
        self.clientdb = ClientDatabase(Session)
        self.dsage_server = DSageServer(self.jobdb,self.workerdb,
                                        self.clientdb, log_level=5)
        self.realm = Realm(self.dsage_server)
        self.p = _SSHKeyPortalRoot(portal.Portal(self.realm))
        self.p.portal.registerChecker(
                            PublicKeyCredentialsCheckerDB(self.clientdb))
        self.client_factory = pb.PBServerFactory(self.p)
        self.hostname = 'localhost'
        self.server = reactor.listenTCP(0, self.client_factory)
        self.port = self.server.getHost().port

        # public key authentication information
        self.username = 'unit_test'
        self.pubkey_file = TEST_PUB_KEY
        self.privkey_file = TEST_PRIV_KEY
        self.data = RANDOM_DATA
        self.public_key = keys.Key.fromString(data=TEST_PUB_KEY)
        self.private_key = keys.Key.fromString(data=TEST_PRIV_KEY)
        self.algorithm = 'rsa'
        self.blob = self.public_key.blob()
        self.signature = self.private_key.sign(self.data)
        self.creds = credentials.SSHPrivateKey(self.username, self.algorithm,
                                               self.blob, self.data,
                                               self.signature)
        pubkey_str = self.public_key.toString('openssh')
        try:
            self.clientdb.del_client(self.username)
            self.clientdb.add_client(self.username, pubkey_str)
        except:
            self.clientdb.add_client(self.username, pubkey_str)

        self.factory = ClientFactory(self._login, (), {})
        self.factory.continueTrying = False
        self.connection = reactor.connectTCP(self.hostname,
                                             self.port,
                                             self.factory)

    def tearDown(self):
        self.jobdb._shutdown()
        os.remove(self.test_db.name)
        self.server.stopListening()
        self.connection.disconnect()


class ClientRemoteCallsTest(RemoteTests):
    """
    Tests of remote procedure calls go here.

    """

    def _catch_failure(self, failure, *args):
        log.msg("Error: ", failure.getErrorMessage())
        log.msg("Traceback: ", failure.printTraceback())

    def testremoteSubmitJob(self):
        """tests perspective_submit_job"""
        jobs = self.create_jobs(1)
        d = self.factory.login(self.creds, None)
        d.addCallback(self._LoginConnected2, jobs)

        return d

    def _LoginConnected2(self, remoteobj, jobs):
        job = jobs[0]
        job.code = "2+2"
        d = remoteobj.callRemote('submit_job', job._reduce())
        d.addCallback(self._got_job_id)

        return d

    def _got_job_id(self, job_id):
        self.assertEquals(type(job_id), str)

    def testremoteSubmitBadJob(self):
        """tests perspective_submit_job"""

        d = self.factory.login(self.creds, None)
        d.addCallback(self._LoginConnected3)

        return d

    def _LoginConnected3(self, remoteobj):
        d = remoteobj.callRemote('submit_job', None)
        d.addErrback(self._gotNoJobID)
        return d

    def _gotNoJobID(self, failure):
        self.assertEquals(BadJobError, failure.check(BadJobError))

    def create_jobs(self, n):
        """This method creates n jobs. """

        jobs = []
        for i in range(n):
            jobs.append(Job(name='unittest', username=self.username))

        return jobs

class WorkerRemoteCallsTest(RemoteTests):
    """
    Tests remote calls for monitors.

    """

    def testremote_get_job(self):
        job = Job()
        job.code = "2+2"
        self.dsage_server.submit_job(job._reduce())
        d = self.factory.login(self.creds, (pb.Referenceable(), hf))
        d.addCallback(self._logged_in)
        d.addCallback(self._get_job)

        return d

    def _logged_in(self, remoteobj):
        self.assert_(remoteobj is not None)

        return remoteobj

    def _get_job(self, remoteobj):
        d = remoteobj.callRemote('get_job')
        d.addCallback(self._got_job)

        return d

    def _got_job(self, jdict):
        self.assertEquals(type(jdict), dict)

    def testremote_job_done(self):
        job = Job()
        job.code = "2+2"
        job_id = self.dsage_server.submit_job(job._reduce())
        jdict = self.dsage_server.get_job()

        d = self.factory.login(self.creds, (pb.Referenceable(), hf))
        d.addCallback(self._logged_in)
        d.addCallback(self._job_done, jdict)

        return d

    def _job_done(self, remoteobj, jdict):
        import time
        job_id = jdict['job_id']
        result = jdict['result']
        d = remoteobj.callRemote('job_done', job_id,
                                 'Nothing.', result, True,
                                 time.time() - time.time())
        d.addCallback(self._done_job)

        return d

    def _done_job(self, job_id):
        self.assertEquals(type(job_id), str)
        jdict = self.dsage_server.get_job_by_id(job_id)
        self.assertEquals(jdict['status'], 'completed')
        self.assertEquals(jdict['output'], 'Nothing.')

    def testremote_job_failed(self):
        job = Job()
        job.code = "2+2"
        job_id = self.dsage_server.submit_job(job._reduce())
        jdict = self.dsage_server.get_job_by_id(job_id)

        d = self.factory.login(self.creds, (pb.Referenceable(), hf))
        d.addCallback(self._logged_in)
        d.addCallback(self._job_failed, jdict)

        return d

    def _job_failed(self, remoteobj, jdict):
        d = remoteobj.callRemote('job_failed', jdict['job_id'], 'Failure')
        d.addCallback(self._failed_job)

        return d

    def _failed_job(self, job_id):
        self.assertEquals(type(job_id), str)
        jdict = self.dsage_server.get_job_by_id(job_id)
        self.assertEquals(jdict['failures'], 1)
        self.assertEquals(jdict['output'], 'Failure')

    def testget_killed_jobs_list(self):
        job = Job()
        job.code = "2+2"
        job.killed = True
        job_id = self.dsage_server.submit_job(job._reduce())
        jdict = self.dsage_server.get_job_by_id(job_id)

        d = self.factory.login(self.creds, (pb.Referenceable(), hf))
        d.addCallback(self._logged_in)
        d.addCallback(self._get_killed_jobs_list)
        d.addCallback(self._got_killed_jobs_list, jdict)

        return d

    def _get_killed_jobs_list(self, remoteobj):
        d = remoteobj.callRemote('get_killed_jobs_list')

        return d

    def _got_killed_jobs_list(self, killed_jobs_list, jdict):
        self.assertEquals(len(killed_jobs_list), 1)
        self.assertEquals(killed_jobs_list[0]['job_id'], jdict['job_id'])

if __name__ == 'main':
    unittest.main()
