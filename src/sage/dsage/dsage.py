"""nodoctests
Distributed SAGE

AUTHORS:
    Yi Qiang (yqiang@gmail.com)
"""

import os

class DistributedSage(object):
    """
    DistributedSage allows you to do distributed computing in SAGE.

    To get up and running quickly, run dsage.setup() to run the
    configuration utility.
    Note that configuration files will be stored in $DOT_SAGE/dsage

    There are three distinct parts of Distributed SAGE:
        Server
            Launch the server with dsage.server()
        Worker
            Launch the worker with dsage.worker()
        Client
            Create the DSage object like this:
                d = DSage()

    Examples:
        This starts a server instance on localhost

        sage: dsage.server()
        This is currently blocking

        Open another sage instance and type:
        sage: dsage.worker()

        This starts a worker connecting the localhost

        Open yet another terminal and type:
        sage: D = DSage()

        This creates a connection to the remote server.  To do a simple
        evaluation, type:

        sage: job1 = D('2+2')

        This sends the job 'print 2+2' to a worker and you can view the
        result by typing:

        sage: print job1

        This is the most basic way of interacting with dsage. To do more
        complicated tasks, you should look at the DistributedFunction class.
        For example, to do distributed integer factorization with ECM, type
        this:

        sage: f = DistributedFactor(P, number, name='my_factor')
        sage: f.start()

        To check the result, do

        sage: print f.result

        To check if it is done, do

        sage: print f.done

    Customization:
        To customize how the worker, server, or client behaves, you can look
        for their respective conf files in DOT_SAGE/dsage.  The configuration
        file should be self explanatory.

    TODO:

    """
    def __init__(self):
        pass

    def start_all(self):
        self.server(blocking=False)
        self.worker(blocking=False)
        from sage.dsage.interface.dsage_interface import BlockingDSage as DSage
        return DSage()

    def server(self, blocking=True):
        """
        This is the server of Distributed SAGE

        Doing dsage.server() will spawn a server process which listens by
        default on ports 8081 and 8082.

        """

        cmd = 'dsage_server.py'
        if not blocking:
            cmd += '&'
        os.system(cmd)

    def worker(self, server=None, port=None, blocking=True):
        """
        This is the worker of Distributed SAGE

        Typing sage.worker() will launch a worker which by default connects to
        localhost on port 8082 to fetch jobs.

        Parameters:
        hostname -- the server you want to connect to
        port -- the port that the server listens on for workers.

        """

        cmd = 'dsage_worker.py'
        cmd += ' %s' % server
        cmd += ' %s' % port

        if not blocking:
            cmd += '&'

        os.system(cmd)

    # This is completely outdated now, only kept for historical reference
    # def console(self):
    #     """
    #     This is the IPython console that allows you submit and view jobs.
    #
    #     Simply type dsage.console() to launch it.  It is a special ipython
    #     console because it has a twisted thread running in the background.
    #     """
    #     # this is overwritten below.
    #     cmd = 'dsage_console.py'
    #     os.system(cmd)

    def setup(self):
        """
        This is the setup utility which helps you configure dsage.

        Type dsage.setup() to run the configuration for the server, worker and
        client.  Alternatively, if you want to run the configuration for just
        one parts, you can launch dsage.setup_server(), dsage.setup_worker() or
        dsage.setup() client

        """

        cmd = 'dsage_setup.py'
        os.system(cmd)

    def setup_server(self):
        """
        This method runs the configuration utility for the server.

        """

        cmd = 'dsage_setup.py server'
        os.system(cmd)

    def setup_worker(self):
        """
        This method runs the configuration utility for the worker.

        """

        cmd = 'dsage_setup.py worker'
        os.system(cmd)

    def setup_client(self):
        """
        This method runs the configuration utility for the client.

        """

        cmd = 'dsage_setup.py client'
        os.system(cmd)

dsage = DistributedSage()
