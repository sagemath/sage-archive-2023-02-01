import unittest
import uuid
import os
import tempfile

from sage.dsage.database.workerdb import WorkerDatabaseSA
from sage.dsage.database.db_config import init_db_sa as init_db
from sage.dsage.misc.hostinfo import HostInfo

class WorkerDatabaseSATestCase(unittest.TestCase):
    test_db = tempfile.NamedTemporaryFile()

    def setUp(self):
        Session = init_db(self.test_db.name)
        self.workerdb = WorkerDatabaseSA(Session)
        self.host_info = HostInfo().host_info
        self.uuid = str(uuid.uuid1())
        self.host_info['uuid'] = self.uuid
        self.host_info['workers'] = 2
        self.workerdb.add_worker(self.host_info)

    def tearDown(self):
        from sqlalchemy.orm import clear_mappers
        self.workerdb.sess.close()
        clear_mappers()
        os.remove(self.test_db.name)

    def testget_worker(self):
        worker = self.workerdb.get_worker(self.uuid)
        self.assertEquals(worker.uuid, self.uuid)

    def testinitial_sate(self):
        worker = self.workerdb.get_worker(self.uuid)
        self.assertEquals(worker.connected, False)

    def testget_worker_count(self):
        self.assertEquals(2, self.workerdb.get_worker_count(False, False))
        self.assertEquals(0, self.workerdb.get_worker_count(True, True))

    def testset_connected(self):
        self.workerdb.set_connected(self.uuid, True)
        worker = self.workerdb.get_worker(self.uuid)
        self.assertEquals(worker.connected, True)