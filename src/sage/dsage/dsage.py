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

dsage = DistributedSage()

