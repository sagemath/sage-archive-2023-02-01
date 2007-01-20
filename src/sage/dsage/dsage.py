import os
class DistributedSage(object):
    def __init__(self):
        pass

    def server(self):
        cmd = 'dsage_server.py'
        os.system(cmd)

    def worker(self):
        cmd = 'dsage_worker.py'
        os.system(cmd)

    def console(self):
        cmd = 'dsage_console.py'
        os.system(cmd)

    def setup():
        cmd = 'dsage_setup.py'
        os.system(cmd)

    def setup_server():
        cmd = 'dsage_setup.py server'
        os.system(cmd)

    def setup_worker():
        cmd = 'dsage_setup.py worker'
        os.system(cmd)

    def setup_client():
        cmd = 'dsage_setup.py client'
        os.system(cmd)

dsage = DistributedSage()
