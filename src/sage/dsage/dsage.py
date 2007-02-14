"""nodoctests
Distributed SAGE

AUTHORS:
    Yi Qiang (yqiang@gmail.com)
"""

import os

class DistributedSage(object):
    r"""
    DistributedSage allows you to do distributed computing in SAGE.

    To get up and running quickly, run dsage.setup() to run the configuration
    utility.
    Note that configuration files will be stored in DOT_SAGE/dsage

    There are three distinct parts of Distributed SAGE:
        Server
            Launch the server with dsage.server()
        Worker
            Launch the worker with dsage.worker()
        Client
            Launch the client with dsage.client()

    Examples:
        This starts a server instance on localhost

        sage: dsage.server()
        This is currently blocking

        Open another sage instance and type

        sage: dsage.worker()

        This starts a worker connecting the localhost

        Open yet another terminal and type:

        sage: dsage.console()

        This drops you into the dsage ipython console

        To do something simple, type:

        sage: P = DSage()

        This creates a connection to the remote server.  To do a simple
        evaluation, type:

        sage: job1 = P.eval('print 2+2', 'testjob')

        This sends the job 'print 2+2' to a worker and you can view the
        result by typing:

        sage: print job1.result

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

    def server(self, blocking=True):
        r"""
        This is the server of Distributed SAGE

        Doing dsage.server() will spawn a server process which listens by
        default on ports 8081 and 8082.

        """

        cmd = 'dsage_server.py'
        if not blocking:
            cmd += '&'
        os.system(cmd)

    def worker(self, hostname=None, port=None, blocking=True):
        r"""
        This is the worker of Distributed SAGE

        Typing sage.worker() will launch a worker which by default connects to
        localhost on port 8082 to fetch jobs.

        Parameters:
        hostname -- the server you want to connect to
        port -- the port that the server listens on for workers.

        """

        cmd = 'dsage_worker.py'
        if isinstance(hostname, str):
            cmd += ' %s' % hostname
        if isinstance(port, int):
            cmd += ' %s' % port

        if not blocking:
            cmd += '&'

        os.system(cmd)

    def console(self):
        r"""
        This is the IPython console that allows you submit and view jobs.

        Simply type dsage.console() to launch it.  It is a special ipython
        console because it has a twisted thread running in the background.
        """
        # this is overwritten below.
        cmd = 'dsage_console.py'
        os.system(cmd)

    def setup(self):
        r"""
        This is the setup utility which helps you configure dsage.

        Type dsage.setup() to run the configuration for the server, worker and
        client.  Alternatively, if you want to run the configuration for just
        one parts, you can launch dsage.setup_server(), dsage.setup_worker() or
        dsage.setup() client

        """
        cmd = 'dsage_setup.py'
        os.system(cmd)

    def setup_server(self):
        r"""
        This method runs the configuration utility for the server.
        """

        cmd = 'dsage_setup.py server'
        os.system(cmd)

    def setup_worker(self):
        r"""
        This method runs the configuration utility for the worker.
        """

        cmd = 'dsage_setup.py worker'
        os.system(cmd)

    def setup_client(self):
        r"""
        This method runs the configuration utility for the client.
        """

        cmd = 'dsage_setup.py client'
        os.system(cmd)

dsage = DistributedSage()

# we have to do it this way, so the proper globals
# get passed to start_dsage_console.
import scripts.dsage_activate
dsage.console = scripts.dsage_activate.start_dsage_console
