import unittest

from sage.dsage.misc.hostinfo import HostInfo

class WorkerTestCase(unittest.TestCase):
    def setUp(self):
        self.worker = Worker(host_info)