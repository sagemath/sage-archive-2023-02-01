r"""
Distributed Sage

Distributed Sage \code{dsage} is a distributed computing framework suitable
for coarse distributed compuatations.

"""
##############################################################################
#
#  DSAGE: Distributed SAGE
#
#       Copyright (C) 2006, 2007 Yi Qiang <yqiang@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#
##############################################################################
import os
import subprocess
from getpass import getuser
import time

import sage.interfaces.cleaner
from sage.misc.all import SAGE_ROOT
from sage.dsage.misc.constants import (DSAGE_DIR, SERVER_LOG, WORKER_LOG,
                                       SERVER_TAC, DSAGE_DB)
from sage.dsage.misc.config import check_dsage_dir
from sage.dsage.misc.misc import find_open_port
import sage.plot.plot

def spawn(cmd, verbose=True, stdout=None, stdin=None):
    """
    Spawns a process and registers it with the SAGE.
    """

    null = open('/dev/null', 'a')
    if stdout is None:
        stdout = null
    if stdin is None:
        stdin = null
    cmdl = cmd.split(' ')
    process = subprocess.Popen(cmdl, shell=False, stdout=stdout, stdin=null)
    sage.interfaces.cleaner.cleaner(process.pid, cmd)
    if verbose:
        print 'Spawned %s (pid = %s)\n' % (' '.join(cmdl), process.pid)

    return process


class DistributedSage(object):
    r"""
    Distributed SAGE allows you to do distributed computing in SAGE.

    To get up and running quickly, run dsage.setup() to run the
    configuration utility.

    Note that configuration files will be stored in the
    directory \code{\$DOT\_SAGE/dsage}.

    QUICK-START

    1.  Launch sage
    2.  Run

        \code{sage: dsage.setup()}

        For a really quick start, just hit ENTER on all questions.
        This will create all the necessary supporting files to get **DSAGE**
        running. It will create the databases, set up a private/public key for
        authentication and create a SSL certificate for the server.
    3.  Launch a server, monitor and get a connection to the server:

        \code{sage: D = dsage.start_all()}

        This will start 2 workers by default. You can change it by passing in
        the ``workers=N`` argument where ``N`` is the number of workers you
        want.
    4.  To do a computation, use D just like any other SAGE interface. For
        example:

        \code{sage: j = D('2+2')}
        \code{sage: j.wait()}
        \code{sage: j}
        \code{4}
    """

    def start_all(self, port=None, workers=2, log_level=0, poll=1.0,
                  authenticate=False, failure_threshold=3,
                  verbose=True, testing=False):
        """
        Start the server and worker and returns a connection to the server.

        """

        from sage.dsage.interface.dsage_interface import BlockingDSage
        from sage.dsage.misc.misc import find_open_port

        if port is None:
            port = find_open_port().next()

        if testing or sage.plot.plot.DOCTEST_MODE:
            testing = True
            self.server(port=port,
                        log_level=5,
                        blocking=False,
                        ssl=False,
                        failure_threshold=failure_threshold,
                        verbose=False,
                        testing=testing)
            self.worker(port=port,
                        workers=workers,
                        log_level=5,
                        ssl=False,
                        blocking=False,
                        poll=0.1,
                        authenticate=authenticate,
                        verbose=False)
        else:
            self.server(port=port,
                        log_level=log_level,
                        blocking=False,
                        failure_threshold=failure_threshold,
                        verbose=verbose)

            self.worker(port=port,
                        workers=workers,
                        log_level=log_level,
                        blocking=False,
                        poll=poll,
                        authenticate=authenticate,
                        verbose=verbose)

        # We want to establish a connection to the server
        tries = 10
        while(tries > 0):
            try:
                import socket
                s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                s.connect(('localhost', port))
                s.close()
                success = True
                break
            except Exception, msg:
                success = False
                tries -= 1
                time.sleep(0.25)
        if not success:
            print 'Could not connect to the server.'
            print 'Error msg from last attempt: %s' % (msg)
            return

        if testing or sage.plot.plot.DOCTEST_MODE:
            d = BlockingDSage(server='localhost', port=port, testing=testing,
                              ssl=False)
        else:
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
            os.kill(self.worker_proc.pid, 9)
            self.worker_proc.wait()
            del self.worker_proc
        except OSError, msg:
            print 'Error killing worker: %s' % msg

    def kill_server(self):
        try:
            os.kill(self.server_proc.pid, 9)
            self.server_proc.wait()
            del self.server_proc
        except OSError, msg:
            print 'Error killing server: %s' % msg

    def server(self, blocking=True, port=None, log_level=0, ssl=True,
               db_file=DSAGE_DB,
               log_file=SERVER_LOG,
               privkey=os.path.join(DSAGE_DIR, 'cacert.pem'),
               cert=os.path.join(DSAGE_DIR, 'pubcert.pem'),
               authenticated_logins=False, failure_threshold=3,
               verbose=True, testing=False, profile=False):
        r"""
        Run the Distributed SAGE server.

        Doing \code{dsage.server()} will spawn a server process which
        listens by default on port 8081.
        """
        open_ports = find_open_port()
        check_dsage_dir()
        cwd = os.getcwd()
        pid_file = 'server.pid'

        def write_tac(tac):
            os.chdir(DSAGE_DIR)
            f = open('dsage_server.tac', 'w')
            f.writelines(tac)
            f.close()

        if testing or sage.plot.plot.DOCTEST_MODE:
            testing = True
            print 'Going into testing mode...'
            # remove the old db to start from a clean slate
            try:
                os.chdir(DSAGE_DIR + '/db')
                os.remove('testing.db')
            except:
                pass
            db_file = 'db/testing.db'

        if port != None:
            server_port = port
        else:
            server_port = open_ports.next()
            open_ports.next()

        tac = SERVER_TAC % (db_file, failure_threshold, ssl, log_level,
                            log_file, privkey, cert, server_port, testing)
        write_tac(tac)

        cmd = 'twistd -d %s --pidfile=%s ' % (DSAGE_DIR, pid_file)
        if profile:
            if verbose:
                print 'Launched with profiling enabled...'
            cmd += '--nothotshot --profile=dsage_server.profile --savestats '
        if blocking:
            cmd += '--nodaemon -y dsage_server.tac'
            cmd += ' | tee -a %s' % (log_file)
            os.system(cmd)
        else:
            try:
                os.remove(pid_file)
            except:
                pass
            cmd += '--logfile=%s -y dsage_server.tac' % (log_file)
            self.server_proc = spawn(cmd, verbose=verbose)
            # Need the following hack since subprocess.Popen reports the wrong
            # pid when launching an application with twistd
            while True:
                try:
                    pid = int(open(pid_file).read())
                    sage.interfaces.cleaner.cleaner(pid, cmd)
                    break
                except:
                    time.sleep(0.1)
                    continue
        os.chdir(cwd)


    def worker(self, server='localhost', port=8081, workers=2, poll=1.0,
               username=getuser(), blocking=True, ssl=True, log_level=0,
               authenticate=True, priority=20,
               privkey=os.path.join(DSAGE_DIR, 'dsage_key'),
               pubkey=os.path.join(DSAGE_DIR, 'dsage_key.pub'),
               log_file=WORKER_LOG,
               verbose=True):
        r"""
        Run the Distributed SAGE worker.

        Typing \code{sage.worker()} will launch a worker which by
        default connects to localhost on port 8081 to fetch jobs.
        """

        check_dsage_dir()
        cmd = ('dsage_worker.py -s %s -p %s -u %s -w %s --poll %s -l %s -f %s '
               + '--privkey=%s --pubkey=%s --priority=%s')
        cmd = cmd % (server, port, username, workers, poll, log_level,
                     log_file, privkey, pubkey, priority)
        if ssl:
            cmd += ' --ssl'
        if authenticate:
            cmd += ' -a'
        if not blocking:
            cmd += ' --noblock'
            cmd = 'python ' + SAGE_ROOT + '/local/bin/' + cmd
            self.worker_proc = spawn(cmd, verbose=verbose)
        else:
            cmd = 'python ' + SAGE_ROOT + '/local/bin/' + cmd
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
