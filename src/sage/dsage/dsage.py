"""nodoctests
Distributed SAGE

AUTHORS:
    Yi Qiang (yqiang@gmail.com)
"""

import os
import subprocess
from getpass import getuser

import sage.interfaces.cleaner
from sage.misc.all import SAGE_ROOT
from sage.dsage.misc.constants import DSAGE_DIR

def spawn(cmd, verbose=True):
    """
    Spawns a process and registers it with the SAGE.
    """

    null = open('/dev/null', 'a')
    cmdl = cmd.split(' ')
    exe = SAGE_ROOT + '/local/bin/' + cmdl[0] # path to the .py file
    cmdl = cmdl[1:]
    proc = [exe] + cmdl + ['&']
    process = subprocess.Popen(proc, shell=False, stdout=null, stdin=null)
    sage.interfaces.cleaner.cleaner(process.pid, cmd)
    if verbose:
        print 'Spawned %s (pid = %s)\n' % (cmd, process.pid)

    return process.pid

class DistributedSage(object):
    r"""
    Distributed SAGE allows you to do distributed computing in SAGE.

    To get up and running quickly, run dsage.setup() to run the
    configuration utility.

    Note that configuration files will be stored in the
    directory \code{\$DOT\_SAGE/dsage}.

    There are three distinct parts of Distributed SAGE:
        Server
            Launch the server with dsage.server()

        Worker
            Launch the worker with dsage.worker()

        Client
            Create the DSage object like this:
                d = DSage()

    EXAMPLES:
    This starts a server instance on localhost:

        sage: dsage.server()

    The dsage server is currently blocking by default.

    Open another sage instance and type:

        sage: dsage.worker()

    This starts a worker connecting the localhost. By default the worker will
    connect to localhost and the port the last server started is listening on.
    All of these settings are configurable via changing
    \code{\$DOT\_SAGE/dsage/worker.conf}

    Open yet another terminal and type:
        sage: D = DSage()

    This creates a connection to the remote server.  To do a simple
    evaluation, type:
        sage: job1 = D('2+2')

    This sends the job '2+2' to a worker and you can view the
    result by typing:

        sage: print job1

    This is the most basic way of interacting with dsage. To do more
    complicated tasks, you should look at the DistributedFunction
    class.  For example, to do distributed integer factorization with
    ECM, type this:

        sage: f = DistributedFactor(P, number, name='my_factor')
        sage: f.start()

    To check the result, do

        sage: print f.result

    To check if it is done, do

        sage: print f.done

    Customization:

        See the DOT\_SAGE/dsage directory.

    """

    def __init__(self):
        pass

    def start_all(self, port=None, workers=2, log_level=0, poll=1.0,
                  anonymous_workers=False, job_failure_threshold=3,
                  verbose=True):
        """
        Start the server and worker and returns a connection to the server.

        """

        from sage.dsage.interface.dsage_interface import BlockingDSage
        from sage.dsage.misc.misc import find_open_port

        if port is None:
            port = find_open_port()
        self.server(port=port, log_level=log_level, blocking=False,
                    job_failure_threshold=job_failure_threshold,
                    verbose=verbose)
        self.worker(port=port, workers=workers, log_level=log_level,
                    blocking=False, poll=poll, anonymous=anonymous_workers,
                    verbose=verbose)

        # We want to establish a connection to the server
        while(True):
            try:
                import socket
                s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                s.connect(('localhost', port))
                s.close()
                break
            except:
                import time
                time.sleep(0.5)

        d = BlockingDSage(server='localhost', port=port)

        return d

    def kill_all(self):
        """
        Kills the server and worker.

        """

        self.kill_worker()
        self.kill_server()

    def kill_worker(self):
        try:
            os.kill(self.worker_pid, 9)
        except OSError, msg:
            print 'Error killing worker: %s' % msg

    def kill_server(self):
        try:
            os.kill(self.server_pid, 9)
        except OSError, msg:
            print 'Error killing server: %s' % msg

    def server(self, blocking=True, port=8081, log_level=0, ssl=True,
               db_file=os.path.join(DSAGE_DIR, 'db', 'dsage.db'),
               log_file=os.path.join(DSAGE_DIR, 'server.log'),
               privkey=os.path.join(DSAGE_DIR, 'cacert.pem'),
               cert=os.path.join(DSAGE_DIR, 'pubcert.pem'),
               stats_file=os.path.join(DSAGE_DIR, 'dsage.xml'),
               anonymous_logins=False,
               job_failure_threshold=3,
               verbose=True):
        r"""
        Run the Distributed SAGE server.

        Doing \code{dsage.server()} will spawn a server process which
        listens by default on port 8081.

        INPUT:
            blocking -- boolean (default: True) -- if False the dsage
                        server will run and you'll still be able to
                        enter commands at the command prompt (though
                        logging will make this hard).
            logfile  -- only used if blocking=True; the default is
                        to log to \file{\$DOT_SAGE/dsage/server.log}

        """

        cmd = 'dsage_server.py -d %s -p %s -l %s -f %s ' + \
                              '-c %s -k %s --jobfailures %s --statsfile=%s'
        cmd = cmd % (db_file, port, log_level, log_file, cert, privkey,
                     job_failure_threshold, stats_file)
        if ssl:
            cmd += ' --ssl'
        if not blocking:
            cmd += ' --noblock'
            self.server_pid = spawn(cmd, verbose=verbose)
        else:
            os.system(cmd)

    def worker(self, server='localhost', port=8081, workers=2, poll=1.0,
               username=getuser(), blocking=True, ssl=True, log_level=0,
               anonymous=False, priority=20,
               privkey=os.path.join(DSAGE_DIR, 'dsage_key'),
               pubkey=os.path.join(DSAGE_DIR, 'dsage_key.pub'),
               log_file=os.path.join(DSAGE_DIR, 'worker.log'),
               verbose=True):
        r"""
        Run the Distributed SAGE worker.

        Typing \code{sage.worker()} will launch a worker which by
        default connects to localhost on port 8081 to fetch jobs.

        INPUT:
            server -- (string, default: None) the server you want to
                      connect to if None, connects to the server
                      specified in .sage/dsage/worker.conf
            port -- (integer, default: None) the port that the server
                      listens on for workers.
            workers -- number of workers to start
            poll -- rate (in seconds) at which the worker pings the server to
                    check for new jobs, this value will increase if the server
                    has no jobs
            username -- username to use
            blocking -- (bool, default: True) whether or not to make a
                        blocking connection.
            ssl -- (bool, default: True) whether or not to use SSL
            log_level -- (int, default: 0) int from 0-5, 5 being most verbose
            anonymous -- (bool, default: False) connect anonymously
            priority -- (int, default: 20) priority of workers
            privkey -- private key
            pubkey -- public key
            log_file -- only used if blocking=True; the default is
                       to log to \file{\$DOT_SAGE/dsage/worker.log}
            verbose -- be more verbose about launching the workers

        """

        cmd = 'dsage_worker.py -s %s -p %s -u %s -w %s --poll %s -l %s -f %s ' + \
                               '--privkey=%s --pubkey=%s --priority=%s '
        cmd = cmd % (server, port, username, workers, poll, log_level,
                     log_file, privkey, pubkey, priority)

        if ssl:
            cmd += ' --ssl'
        if anonymous:
            cmd += ' -a'
        if not blocking:
            cmd += ' --noblock'
            self.worker_pid = spawn(cmd, verbose=verbose)
        else:
            os.system(cmd)

    def setup(self, template=None):
        r"""
        This is the setup utility which helps you configure dsage.

        Type \code{dsage.setup()} to run the configuration for the server,
        worker and client.  Alternatively, if you want to run the
        configuration for just one parts, you can launch
        \code{dsage.setup_server()}, \code{dsage.setup\_worker()}
        or \code{dsage.setup()}.

        """

        from sage.dsage.scripts.dsage_setup import setup
        setup(template=template)

    def setup_server(self, *args):
        """
        This method runs the configuration utility for the server.

        """

        from sage.dsage.scripts.dsage_setup import setup_server
        setup_server(*args)

    def setup_worker(self):
        """
        This method runs the configuration utility for the worker.

        """

        from sage.dsage.scripts.dsage_setup import setup_worker
        setup_worker()

    def setup_client(self):
        """
        This method runs the configuration utility for the client.

        """

        from sage.dsage.scripts.dsage_setup import setup_client
        setup_client()

dsage = DistributedSage()
