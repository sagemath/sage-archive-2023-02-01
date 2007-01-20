import os
class DistributedSage(object):
    def __init__(self):
        pass

    def server(self, blocking=False):
        cmd = 'dsage_server.py'
        if not blocking:
            cmd += '&'
        os.system(cmd)

    def worker(self, blocking=False):
        cmd = 'dsage_worker.py'
        if not blocking:
            cmd += '&'
        os.system(cmd)

    def console(self):
        cmd = 'dsage_console.py'
        os.system(cmd)

    def setup(self):
        cmd = 'dsage_setup.py'
        os.system(cmd)

    def setup_server(self):
        cmd = 'dsage_setup.py server'
        os.system(cmd)

    def setup_worker(self):
        cmd = 'dsage_setup.py worker'
        os.system(cmd)

    def setup_client(self):
        cmd = 'dsage_setup.py client'
        os.system(cmd)

dsage = DistributedSage()
