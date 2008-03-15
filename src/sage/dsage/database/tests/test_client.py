import unittest

from sage.dsage.database.client import Client

class ClientTestCase(unittest.TestCase):
    username = 'tester'
    pub_key = 'tester_key'
    def setUp(self):
        self.c = Client(self.username, self.pub_key)

    def testInitClient(self):
        self.assertEquals(self.c.username, self.username)
        self.assertEquals(self.c.public_key, self.pub_key)

    def testGetUsername(self):
        self.assertEquals(self.c.get_username(), self.username)

    def testGetPublicKey(self):
        self.assertEquals(self.c.get_public_key(), self.pub_key)

    def testSetConnected(self):
        self.c.connected = False
        self.assertEquals(self.c.is_connected(), False)
        self.c.connected = True
        self.assertEquals(self.c.is_connected(), True)

    def testSetEnabled(self):
        self.c.enabled = False
        self.assertEquals(self.c.is_enabled(), False)
        self.c.enabled = True
        self.assertEquals(self.c.is_enabled(), True)
