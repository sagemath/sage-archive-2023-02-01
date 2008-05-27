#!/usr/bin/env python
############################################################################
#
#   DSAGE: Distributed SAGE
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
############################################################################
__docformat__ = "restructuredtext en"

import sys
import os
import cPickle
import zlib
import pexpect
import datetime
from math import ceil
from getpass import getuser

from twisted.spread import pb
from twisted.internet import reactor, defer, error, task
from twisted.python import log
from twisted.spread import banana
banana.SIZE_LIMIT = 100*1024*1024 # 100 MegaBytes

from gnutls.constants import *
from gnutls.crypto import *
from gnutls.errors import *
from gnutls.interfaces.twisted import X509Credentials

from sage.interfaces.sage0 import Sage
from sage.misc.preparser import preparse_file

from sage.dsage.database.job import Job, expand_job
from sage.dsage.misc.hostinfo import HostInfo
from sage.dsage.errors.exceptions import NoJobException
from sage.dsage.twisted.pb import ClientFactory
from sage.dsage.misc.constants import DELIMITER
from sage.dsage.misc.constants import DSAGE_DIR
from sage.dsage.misc.constants import TMP_WORKER_FILES
from sage.dsage.misc.misc import random_str, get_uuid

START_MARKER = '\x01r\x01e'
END_MARKER = '\x01r\x01b'
LOG_PREFIX = "[Worker %s] "

class Worker(object):
    """
    Workers perform the computation of dsage jobs.

    """

    def __init__(self, remoteobj, id, log_level=0, poll=1.0):
        """
        :type remoteobj: remoteobj
        :param remoteobj: Reference to the remote dsage server

        :type id: integer
        :param id: numerical identifier of worker

        :type log_level: integer
        :param log_level: log level, higher means more verbose

        :type poll: integer
        :param poll: rate (in seconds) a worker talks to the server

        """

        self.remoteobj = remoteobj
        self.id = id
        self.free = True
        self.job = None
        self.log_level = log_level
        self.poll_rate = poll
        self.checker_task = task.LoopingCall(self.check_work)
        self.checker_timeout = 0.5
        self.got_output = False
        self.job_start_time = None
        self.orig_poll = poll
        self.start()

    def _catch_failure(self, failure):
        log.msg("Error: ", failure.getErrorMessage())
        log.msg("Traceback: ", failure.printTraceback())

    def _increase_poll_rate(self):
        if self.poll_rate >= 15: # Cap the polling interval to 15 seconds
            self.poll_rate = 15
            if self.log_level > 3:
                log.msg('[Worker %s] Capping poll rate to %s'
                         % (self.id, self.poll_rate))
        else:
            self.poll_rate = ceil(self.poll_rate * 1.5)
            if self.log_level > 3:
                log.msg('[Worker %s] Increased polling rate to %s'
                        % (self.id, self.poll_rate))

    def get_job(self):
        try:
            if self.log_level > 3:
                log.msg(LOG_PREFIX % self.id +  'Getting job...')
            d = self.remoteobj.callRemote('get_job')
        except Exception, msg:
            log.msg(msg)
            log.msg(LOG_PREFIX % self.id +  'Disconnected...')
            self._increase_poll_rate()
            reactor.callLater(self.poll_rate, self.get_job)
            return
        d.addCallback(self.gotJob)
        d.addErrback(self.noJob)

        return d

    def gotJob(self, jdict):
        """
        callback for the remoteobj's get_job method.

        :type jdict: dict
        :param jdict: job dictionary

        """

        if self.log_level > 1:
            if jdict is None:
                log.msg(LOG_PREFIX % self.id + 'No new job.')
        if self.log_level > 3:
            if jdict is not None:
                log.msg(LOG_PREFIX % self.id + 'Got Job: %s' % jdict)
        self.job = expand_job(jdict)
        if not isinstance(self.job, Job):
            raise NoJobException
        try:
            self.poll_rate = self.orig_poll
            self.doJob(self.job)
        except Exception, msg:
            log.msg(msg)
            self.report_failure(msg)
            self.restart()

    def job_done(self, output, result, completed, cpu_time):
        """
        Reports to the server that a job has finished. It also reports partial
        completeness by presenting the server with new output.

        Parameters:
        :type output: string
        :param output: output of command (to sys.stdout)

        :type result: python pickle
        :param result: result of the job

        :type completed: bool
        :param completed: whether or not the job is finished

        :type cpu_time: string
        :param cpu_time: how long the job took

        """

        job_id = self.job.job_id
        wait = 5.0
        try:
            d = self.remoteobj.callRemote('job_done', job_id, output, result,
                                          completed, cpu_time)
        except Exception, msg:
            log.msg('Error trying to submit job status...')
            log.msg('Retrying to submit again in %s seconds...' % wait)
            log.err(msg)
            reactor.callLater(wait, self.job_done, output, result,
                              completed, cpu_time)
            d = defer.Deferred()
            d.errback(error.ConnectionLost())
            return d

        if completed:
            log.msg('[Worker %s] Finished job %s' % (self.id, job_id))
            self.restart()

        return d


    def noJob(self, failure):
        """
        Errback that catches the NoJobException.

        :type failure: twisted.python.failure
        :param failure: a twisted failure object

        """

        if failure.check(NoJobException):
            if self.log_level > 1:
                msg = 'Sleeping for %s seconds' % self.poll_rate
                log.msg(LOG_PREFIX % self.id + msg)
            self._increase_poll_rate()
            reactor.callLater(self.poll_rate, self.get_job)
        else:
            log.msg("Error: ", failure.getErrorMessage())
            log.msg("Traceback: ", failure.printTraceback())

    def setup_tmp_dir(self, job):
        """
        Creates the temporary directory for the worker.

        :type job: sage.dsage.database.job.Job
        :param job: a Job object

        """

        cur_dir = os.getcwd() # keep a reference to the current directory
        tmp_job_dir = os.path.join(TMP_WORKER_FILES, job.job_id)
        if not os.path.isdir(TMP_WORKER_FILES):
            os.mkdir(TMP_WORKER_FILES)
        if not os.path.isdir(tmp_job_dir):
            os.mkdir(tmp_job_dir)
        os.chdir(tmp_job_dir)
        self.sage.eval("os.chdir('%s')" % tmp_job_dir)

        return tmp_job_dir

    def extract_and_load_job_data(self, job):
        """
        Extracts all the data that is in a job object.

        :type job: sage.dsage.database.job.Job
        :param job: a Job object

        """

        if isinstance(job.data, list):
            if self.log_level > 2:
                msg = 'Extracting job data...'
                log.msg(LOG_PREFIX % self.id + msg)
            try:
                for var, data, kind in job.data:
                    try:
                        data = zlib.decompress(data)
                    except Exception, msg:
                        log.msg(msg)
                        continue
                    if kind == 'file':
                        data = preparse_file(data, magic=True, do_time=False,
                                             ignore_prompts=False)
                        f = open(var, 'wb')
                        f.write(data)
                        f.close()
                        if self.log_level > 2:
                            msg = 'Extracted %s' % f
                            log.msg(LOG_PREFIX % self.id + msg)
                        self.sage.eval("execfile('%s')" % var)
                    if kind == 'object':
                        fname = var + '.sobj'
                        if self.log_level > 2:
                            log.msg('Object to be loaded: %s' % fname)
                        f = open(fname, 'wb')
                        f.write(data)
                        f.close()
                        self.sage.eval("%s = load('%s')" % (var, fname))
                        if self.log_level > 2:
                            msg = 'Loaded %s' % fname
                            log.msg(LOG_PREFIX % self.id + msg)
            except Exception, msg:
                log.msg(LOG_PREFIX % self.id + msg)

    def write_job_file(self, job):
        """
        Writes out the job file to be executed to disk.

        :type job: sage.dsage.database.job.Job
        :param job: A Job object

        """

        parsed_file = preparse_file(job.code, magic=True,
                                    do_time=False, ignore_prompts=False)

        job_filename = str(job.name) + '.py'
        job_file = open(job_filename, 'w')
        BEGIN = "print '%s'\n\n" % (START_MARKER)
        END = "print '%s'\n\n" % (END_MARKER)
        GO_TO_TMP_DIR = """os.chdir('%s')\n""" % self.tmp_job_dir
        SAVE_TIME = """save((time.time()-dsage_start_time), 'cpu_time.sobj', compress=False)\n"""
        SAVE_RESULT = """try:
    save(DSAGE_RESULT, 'result.sobj', compress=True)
except:
    save('No DSAGE_RESULT', 'result.sobj', compress=True)
"""
        job_file.write("alarm(%s)\n\n" % (job.timeout))
        job_file.write("import time\n\n")
        job_file.write(BEGIN)
        job_file.write('dsage_start_time = time.time()\n')
        job_file.write(parsed_file)
        job_file.write("\n\n")
        job_file.write(END)
        job_file.write("\n")
        job_file.write(GO_TO_TMP_DIR)
        job_file.write(SAVE_RESULT)
        job_file.write(SAVE_TIME)
        job_file.close()
        if self.log_level > 2:
            log.msg('[Worker: %s] Wrote job file. ' % (self.id))

        return job_filename

    def doJob(self, job):
        """
        Executes a job

        :type job: sage.dsage.database.job.Job
        :param job: A Job object

        """

        log.msg(LOG_PREFIX % self.id + 'Starting job %s ' % job.job_id)

        self.free = False
        self.got_output = False
        d = defer.Deferred()

        try:
            self.checker_task.start(self.checker_timeout, now=False)
        except AssertionError:
            self.checker_task.stop()
            self.checker_task.start(self.checker_timeout, now=False)
        if self.log_level > 2:
            log.msg(LOG_PREFIX % self.id + 'Starting checker task...')

        self.tmp_job_dir = self.setup_tmp_dir(job)
        self.extract_and_load_job_data(job)

        job_filename = self.write_job_file(job)

        f = os.path.join(self.tmp_job_dir, job_filename)
        self.sage._send("execfile('%s')" % (f))
        self.job_start_time = datetime.datetime.now()
        if self.log_level > 2:
            msg = 'File to execute: %s' % f
            log.msg(LOG_PREFIX % self.id + msg)

        d.callback(True)

    def reset_checker(self):
        """
        Resets the output/result checker for the worker.

        """

        if self.checker_task.running:
            self.checker_task.stop()
        self.checker_timeout = 1.0
        self.checker_task = task.LoopingCall(self.check_work)

    def check_work(self):
        """
        check_work periodically polls workers for new output. The period is
        determined by an exponential back off algorithm.

        This figures out whether or not there is anything new output that we
        should submit to the server.

        """

        if self.sage == None:
            return
        if self.job == None or self.free == True:
            if self.checker_task.running:
                self.checker_task.stop()
            return
        if self.log_level > 1:
            msg = 'Checking job %s' % self.job.job_id
            log.msg(LOG_PREFIX % self.id + msg)
        os.chdir(self.tmp_job_dir)
        try:
            # foo, output, new = self.sage._so_far()
            # This sucks and is a very bad way to tell when a calculation is
            # finished
            done, new = self.sage._get()
            # If result.sobj exists, our calculation is done
            result = open('result.sobj', 'rb').read()
            done = True
        except RuntimeError, msg: # Error in calling worker.sage._so_far()
            done = False
            if self.log_level > 1:
                log.msg(LOG_PREFIX % self.id + 'RuntimeError: %s' % msg)
                log.msg("Don't worry, the RuntimeError above " +
                        "is a non-fatal SAGE failure")
            self.increase_checker_task_timeout()
            return
        except IOError, msg: # File does not exist yet
            done = False

        if done:
            try:
                cpu_time = cPickle.loads(open('cpu_time.sobj', 'rb').read())
            except IOError:
                cpu_time = -1
            self.free = True
            self.reset_checker()
        else:
            result = cPickle.dumps('Job not done yet.', 2)
            cpu_time = None

        if self.check_failure(new):
            self.report_failure(new)
            self.restart()
            return

        sanitized_output = self.clean_output(new)
        if self.log_level > 3:
            print 'Output before sanitizing: \n' , sanitized_output
        if self.log_level > 3:
            print 'Output after sanitizing: \n', sanitized_output
        if sanitized_output == '' and not done:
            self.increase_checker_task_timeout()
        else:
            d = self.job_done(sanitized_output, result, done, cpu_time)
            d.addErrback(self._catch_failure)

    def report_failure(self, failure):
        """
        Reports failure of a job.

        :type failure: twisted.python.failure
        :param failure: A twisted failure object

        """

        msg = 'Job %s failed!' % (self.job.job_id)
        import shutil
        failed_dir = self.tmp_job_dir + '_failed'
        if os.path.exists(failed_dir):
            shutil.rmtree(failed_dir)
        shutil.move(self.tmp_job_dir, failed_dir)
        log.msg(LOG_PREFIX % self.id + msg)
        log.msg('Traceback: \n%s' % failure)
        d = self.remoteobj.callRemote('job_failed', self.job.job_id, failure)
        d.addErrback(self._catch_failure)

        return d

    def increase_checker_task_timeout(self):
        """
        Quickly decreases the number of times a worker checks for output

        """

        if self.checker_task.running:
            self.checker_task.stop()

        self.checker_timeout = self.checker_timeout * 1.5
        if self.checker_timeout > 300.0:
            self.checker_timeout = 300.0
        self.checker_task = task.LoopingCall(self.check_work)
        self.checker_task.start(self.checker_timeout, now=False)
        if self.log_level > 0:
            msg = 'Checking output again in %s' % self.checker_timeout
            log.msg(LOG_PREFIX % self.id + msg)

    def clean_output(self, sage_output):
        """
        clean_output attempts to clean up the output string from sage.

        :type sage_output: string
        :param sage_output: sys.stdout output from the child sage instance

        """

        begin = sage_output.find(START_MARKER)
        if begin != -1:
            self.got_output = True
            begin += len(START_MARKER)
        else:
            begin = 0
        end = sage_output.find(END_MARKER)
        if end != -1:
            end -= 1
        else:
            if not self.got_output:
                end = 0
            else:
                end = len(sage_output)
        output = sage_output[begin:end]
        output = output.strip()
        output = output.replace('\r', '')

        if ('execfile' in output or 'load' in output) and self.got_output:
            output = ''

        return output

    def check_failure(self, sage_output):
        """
        Checks for signs of exceptions or errors in the output.

        :type sage_output: string
        :param sage_output: output from the sage instance

        """

        if sage_output == None:
            return False
        else:
            sage_output = ''.join(sage_output)

        if 'Traceback' in sage_output:
            return True
        elif 'Error' in sage_output:
            return True
        else:
            return False

    def kill_sage(self):
        """
        Try to hard kill the SAGE instance.

        """

        try:
            self.sage.quit()
            del self.sage
        except Exception, msg:
            pid = self.sage.pid()
            cmd = 'kill -9 %s' % pid
            os.system(cmd)
            log.msg(msg)

    def stop(self, hard_reset=False):
        """
        Stops the current worker and resets it's internal state.

        :type hard_reset: boolean
        :param hard_reset: Specifies whether to kill -9 the sage instances

        """

        # Set status to free and delete any current jobs we have
        self.free = True
        self.job = None

        if hard_reset:
            log.msg(LOG_PREFIX % self.id + 'Performing hard reset.')
            self.kill_sage()
        else: # try for a soft reset
            INTERRUPT_TRIES = 20
            timeout = 0.3
            e = self.sage._expect
            try:
                for i in range(INTERRUPT_TRIES):
                    self.sage._expect.sendline('q')
                    self.sage._expect.sendline(chr(3))  # send ctrl-c
                    try:
                        e.expect(self.sage._prompt, timeout=timeout)
                        success = True
                        break
                    except (pexpect.TIMEOUT, pexpect.EOF), msg:
                        success = False
                        if self.log_level > 3:
                            msg = 'Interrupting SAGE (try %s)' % i
                            log.msg(LOG_PREFIX % self.id + msg)
            except Exception, msg:
                success = False
                log.msg(msg)
                log.msg(LOG_PREFIX % self.id + "Performing hard reset.")

            if not success:
                self.kill_sage()
            else:
                self.sage.reset()

    def start(self):
        """
        Starts a new worker if it does not exist already.

        """

        log.msg('[Worker %s] Started...' % (self.id))
        if not hasattr(self, 'sage'):
            if self.log_level > 3:
                logfile = DSAGE_DIR + '/%s-pexpect.log' % self.id
                self.sage = Sage(maxread=1, logfile=logfile, python=True)
            else:
                self.sage = Sage(maxread=1, python=True)
            try:
                self.sage._start(block_during_init=True)
            except RuntimeError, msg: # Could not start SAGE
                print msg
                print 'Failed to start a worker, probably Expect issues.'
                reactor.stop()
                sys.exit(-1)
        E = self.sage.expect()
        E.sendline('\n')
        E.expect('>>>')
        cmd = 'from sage.all import *;'
        cmd += 'from sage.all_notebook import *;'
        cmd += 'import sage.server.support as _support_; '
        cmd += 'import time;'
        cmd += 'import os;'
        E.sendline(cmd)

        if os.uname()[0].lower() == 'linux':
            try:
                self.base_mem = int(self.sage.get_memory_usage())
            except:
                pass

        self.get_job()

    def restart(self):
        """
        Restarts the current worker.

        """

        log.msg('[Worker: %s] Restarting...' % (self.id))

        if hasattr(self, 'base_mem'):
            try:
                cur_mem = int(self.sage.get_memory_usage())
            except:
                cur_mem = 0
        try:
            if hasattr(self, 'base_mem'):
                if cur_mem >= (2 * self.base_mem):
                    self.stop(hard_reset=True)
            else:
                from sage.dsage.misc.misc import timedelta_to_seconds
                delta = datetime.datetime.now() - self.job_start_time
                secs = timedelta_to_seconds(delta)
                if secs >= (3*60): # more than 3 minutes, do a hard reset
                    self.stop(hard_reset=True)
                else:
                    self.stop(hard_reset=False)
        except TypeError:
            self.stop(hard_reset=True)
        self.job_start_time = None
        self.start()
        self.reset_checker()


class Monitor(pb.Referenceable):
    """
    Monitors control workers.
    They are able to shutdown workers and spawn them, as well as check on
    their status.

    """

    def __init__(self, server='localhost', port=8081, username=getuser(),
                 ssl=True, workers=2, authenticate=False, priority=20,
                 poll=1.0, log_level=0,
                 log_file=os.path.join(DSAGE_DIR, 'worker.log'),
                 pubkey_file=None, privkey_file=None):
        """
        :type server: string
        :param server: hostname of remote server

        :type port: integer
        :param port: port of remote server

        :type username: string
        :param username: username to use for authentication

        :type ssl: boolean
        :param ssl: specify whether or not to use SSL for the connection

        :type workers: integer
        :param workers: specifies how many workers to launch

        :type authenticate: boolean
        :param authenticate: specifies whether or not to authenticate

        :type priority: integer
        :param priority: specifies the UNIX priority of the workers

        :type poll: float
        :param poll: specifies how fast workers talk to the server in seconds

        :type log_level: integer
        :param log_level: specifies verbosity of logging, higher equals more

        :type log_file: string
        :param log_file: specifies the location of the log_file

        """

        self.server = server
        self.port = port
        self.username = username
        self.ssl = ssl
        self.workers = workers
        self.authenticate = authenticate
        self.priority = priority
        self.poll_rate = poll
        self.log_level = log_level
        self.log_file = log_file
        self.pubkey_file = pubkey_file
        self.privkey_file = privkey_file

        self.remoteobj = None
        self.connected = False
        self.reconnecting = False
        self.worker_pool = None
        self.sleep_time = 1.0

        self.host_info = HostInfo().host_info

        self.host_info['uuid'] = get_uuid()
        self.host_info['workers'] = self.workers
        self.host_info['username'] = self.username

        self._startLogging(self.log_file)

        try:
            os.nice(self.priority)
        except OSError, msg:
            log.msg('Error setting priority: %s' % (self.priority))
            pass
        if self.authenticate:
            from twisted.cred import credentials
            from twisted.conch.ssh import keys
            self.DATA =  random_str(500)
            # public key authentication information
            self.pubkey = keys.Key.fromFile(self.pubkey_file)
            # try getting the private key object without a passphrase first
            try:
                self.privkey = keys.Key.fromFile(self.privkey_file)
            except keys.BadKeyError:
                pphrase = self._getpassphrase()
                self.privkey = keys.Key.fromFile(self.privkey_file,
                                                  passphrase=pphrase)
            self.algorithm = 'rsa'
            self.blob = self.pubkey.blob()
            self.data = self.DATA
            self.signature = self.privkey.sign(self.data)
            self.creds = credentials.SSHPrivateKey(self.username,
                                                   self.algorithm,
                                                   self.blob,
                                                   self.data,
                                                   self.signature)

    def _startLogging(self, log_file):
        """
        :type log_file: string
        :param log_file: file name to log to

        """

        if log_file == 'stdout':
            log.startLogging(sys.stdout)
            log.msg('WARNING: Only loggint to stdout!')
        else:
            worker_log = open(log_file, 'a')
            log.startLogging(sys.stdout)
            log.startLogging(worker_log)
            log.msg("Logging to file: ", log_file)

    def _getpassphrase(self):
        import getpass
        passphrase = getpass.getpass('Passphrase (Hit enter for None): ')

        return passphrase

    def _connected(self, remoteobj):
        """
        Callback for connect.

        :type remoteobj: remote object
        :param remoteobj: remote obj

        """

        self.remoteobj = remoteobj
        self.remoteobj.notifyOnDisconnect(self._disconnected)
        self.connected = True

        if self.worker_pool == None: # Only pool workers the first time
            self.pool_workers(self.remoteobj)
        else:
            for worker in self.worker_pool:
                worker.remoteobj = self.remoteobj # Update workers
                if worker.job == None:
                    worker.restart()

    def _disconnected(self, remoteobj):
        """
        :type remoteobj: remote object
        :param remoteobj: remote obj

        """

        log.msg('Closed connection to the server.')
        self.connected = False

    def _got_killed_jobs(self, killed_jobs):
        """
        Callback for check_killed_jobs.

        :type killed_jobs: dict
        :param killed_jobs: dict of job jdicts which were killed

        """

        if killed_jobs == None:
            return
        killed_jobs = [expand_job(jdict) for jdict in killed_jobs]
        for worker in self.worker_pool:
            if worker.job is None:
                continue
            if worker.free:
                continue
            for job in killed_jobs:
                if job is None or worker.job is None:
                    continue
                if worker.job.job_id == job.job_id:
                    msg = 'Processing killed job, restarting...'
                    log.msg(LOG_PREFIX % worker.id + msg)
                    worker.restart()

    def _retryConnect(self):
        log.msg('[Monitor] Disconnected, reconnecting in %s' % (5.0))
        if not self.connected:
            reactor.callLater(5.0, self.connect)

    def _catchConnectionFailure(self, failure):
        log.msg("Error: ", failure.getErrorMessage())
        log.msg("Traceback: ", failure.printTraceback())
        self._disconnected(None)

    def _catch_failure(self, failure):
        log.msg("Error: ", failure.getErrorMessage())
        log.msg("Traceback: ", failure.printTraceback())

    def connect(self):
        """
        This method connects the monitor to a remote PB server.

        """

        if self.connected: # Don't connect multiple times
            return

        self.factory = ClientFactory(self._login, (), {})
        cred = None
        if self.ssl:
            cred = X509Credentials()
            reactor.connectTLS(self.server, self.port, self.factory, cred)
        else:
            reactor.connectTCP(self.server, self.port, self.factory)

        log.msg(DELIMITER)
        log.msg('DSAGE Worker')
        log.msg('Started with PID: %s' % (os.getpid()))
        log.msg('Connecting to %s:%s' % (self.server, self.port))
        if cred is not None:
            log.msg('Using SSL: True')
        else:
            log.msg('Using SSL: False')
        log.msg(DELIMITER)

    def _login(self, *args, **kwargs):
        if self.authenticate:
            log.msg('Connecting as authenticated worker...\n')
            d = self.factory.login(self.creds, (self, self.host_info))
        else:
            from twisted.cred.credentials import Anonymous
            log.msg('Connecting as unauthenticated worker...\n')
            d = self.factory.login(Anonymous(), (self, self.host_info))
        d.addCallback(self._connected)
        d.addErrback(self._catchConnectionFailure)

        return d

    def pool_workers(self, remoteobj):
        """
        Creates the worker pool.

        """

        log.msg('[Monitor] Starting %s workers...' % (self.workers))
        self.worker_pool = [Worker(remoteobj, x, self.log_level,
                            self.poll_rate)
                            for x in range(self.workers)]


    def remote_set_uuid(self, uuid):
        """
        Sets the workers uuid.
        This is called by the server.

        """

        from sage.dsage.misc.misc import set_uuid
        set_uuid(uuid)


    def remote_calc_score(self, script):
        """
        Calculuates the worker score.

        :type script: string
        :param script: script to score the worker

        """

        from sage.dsage.misc.misc import exec_wrs

        return exec_wrs(script)


    def remote_kill_job(self, job_id):
        """
        Kills the job given the job id.

        :type job_id: string
        :param job_id: the unique job identifier.

        """

        print 'Killing %s' % (job_id)
        for worker in self.worker_pool:
            if worker.job != None:
                if worker.job.job_id == job_id:
                    worker.restart()


def usage():
    """
    Prints usage help.

    """

    from optparse import OptionParser

    usage = ['usage: %prog [options]\n',
              'Bug reports to <yqiang@gmail.com>']
    parser = OptionParser(usage=''.join(usage))
    parser.add_option('-s', '--server',
                      dest='server',
                      default='localhost',
                      help='hostname. Default is localhost')
    parser.add_option('-p', '--port',
                      dest='port',
                      type='int',
                      default=8081,
                      help='port to connect to. default=8081')
    parser.add_option('--poll',
                      dest='poll',
                      type='float',
                      default=5.0,
                      help='poll rate before checking for new job. default=5')
    parser.add_option('-a', '--authenticate',
                      dest='authenticate',
                      default=False,
                      action='store_true',
                      help='Connect as authenticate worker. default=True')
    parser.add_option('-f', '--logfile',
                      dest='logfile',
                      default=os.path.join(DSAGE_DIR, 'worker.log'),
                      help='log file')
    parser.add_option('-l', '--loglevel',
                      dest='loglevel',
                      type='int',
                      default=0,
                      help='log level. default=0')
    parser.add_option('--ssl',
                      dest='ssl',
                      action='store_true',
                      default=False,
                      help='enable or disable ssl')
    parser.add_option('--privkey',
                      dest='privkey_file',
                      default=os.path.join(DSAGE_DIR, 'dsage_key'),
                      help='private key file. default = ' +
                           '~/.sage/dsage/dsage_key')
    parser.add_option('--pubkey',
                      dest='pubkey_file',
                      default=os.path.join(DSAGE_DIR, 'dsage_key.pub'),
                      help='public key file. default = ' +
                           '~/.sage/dsage/dsage_key.pub')
    parser.add_option('-w', '--workers',
                      dest='workers',
                      type='int',
                      default=2,
                      help='number of workers. default=2')
    parser.add_option('--priority',
                      dest='priority',
                      type='int',
                      default=20,
                      help='priority of workers. default=20')
    parser.add_option('-u', '--username',
                      dest='username',
                      default=getuser(),
                      help='username')
    parser.add_option('--noblock',
                      dest='noblock',
                      action='store_true',
                      default=False,
                      help='tells that the server was ' +
                           'started in blocking mode')
    (options, args) = parser.parse_args()

    return options

def main():
    options = usage()
    SSL = options.ssl
    monitor = Monitor(server=options.server, port=options.port,
                      username=options.username, ssl=SSL,
                      workers=options.workers,
                      authenticate=options.authenticate,
                      priority=options.priority, poll=options.poll,
                      log_file=options.logfile,
                      log_level=options.loglevel,
                      pubkey_file=options.pubkey_file,
                      privkey_file=options.privkey_file)
    monitor.connect()
    try:
        if options.noblock:
            reactor.run(installSignalHandlers=0)
        else:
            reactor.run(installSignalHandlers=1)
    except:
        log.msg('Error starting the twisted reactor, exiting...')
        sys.exit()

if __name__ == '__main__':
    usage()
    main()
