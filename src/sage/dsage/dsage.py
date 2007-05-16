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
    Spawns a process and registers it with the SAGE cleaner.

    """

    null = open('/dev/null', 'a')
    proc = '%s/%s' % (SAGE_ROOT + '/local/bin', cmd)
    process = subprocess.Popen(proc, shell=True, stdout=null, stderr=null)
    sage.interfaces.cleaner.cleaner(process.pid, cmd)
    if verbose:
        print 'Spawned %s (pid = %s)\n' % (cmd, process.pid)

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

        To customize how the worker, server, or client behaves, you
        can look for their respective conf files in DOT_SAGE/dsage.
        The configuration file should be self explanatory.

    """

    def __init__(self):
        pass

    def start_all(self, port=8081, log_level=0, verbose=True):
        """
        Start the server and worker and returns a connection to the server.

        """

        from sage.dsage.interface.dsage_interface import BlockingDSage

        self.server(port=port, log_level=log_level, blocking=False,
                    verbose=verbose)
        self.worker(port=port, log_level=log_level, blocking=False,
                    verbose=verbose)

        d = BlockingDSage(server='localhost', port=port)

        return d

    def server(self, blocking=True, port=8081, log_level=0, ssl=True,
               db_file=os.path.join(DSAGE_DIR, 'db', 'dsage.db'),
               log_file=os.path.join(DSAGE_DIR, 'server.log'),
               privkey=os.path.join(DSAGE_DIR, 'cacert.pem'),
               cert=os.path.join(DSAGE_DIR, 'pubcert.pem'),
               stats_file=os.path.join(DSAGE_DIR, 'dsage.xml'),
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
                        to log to $DOT_SAGE/dsage/server.log

        """

        cmd = 'dsage_server.py -d %s -p %s -l %s -f %s ' + \
                              '-c %s -k %s --statsfile=%s --ssl=%s'
        cmd = cmd % (db_file, port, log_level, log_file, cert, privkey,
                     stats_file, ssl)
        if not blocking:
            cmd += ' --noblock'
            spawn(cmd, verbose=verbose)
        else:
            os.system(cmd)

    def worker(self, server='localhost', port=8081, workers=2, delay=5.0,
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
            blocking -- (bool, default: True) whether or not to make a
                        blocking connection.
            logfile -- only used if blocking=True; the default is
                       to log to $DOT_SAGE/dsage/worker.log
        """

        cmd = 'dsage_worker.py -s %s -p %s -u %s -w %s -d %s -l %s -f %s ' + \
                               '--privkey %s --pubkey %s --priority %s'
        cmd = cmd % (server, port, username, workers, delay, log_level,
                     log_file, privkey, pubkey, priority)

        if not blocking:
            cmd += ' --noblock'
            spawn(cmd, verbose=verbose)
        else:
            os.system(cmd)

    def setup(self):
        r"""
        This is the setup utility which helps you configure dsage.

        Type \code{dsage.setup()} to run the configuration for the server,
        worker and client.  Alternatively, if you want to run the
        configuration for just one parts, you can launch
        \code{dsage.setup_server()}, \code{dsage.setup\_worker()}
        or \code{dsage.setup()}.

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
