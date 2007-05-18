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

import sys
import os
import cPickle
import zlib
import pexpect
import datetime
import uuid
from getpass import getuser

from twisted.spread import pb
from twisted.internet import reactor, defer, error, task
from twisted.python import log
from twisted.spread import banana
banana.SIZE_LIMIT = 100*1024*1024 # 100 MegaBytes

from sage.interfaces.sage0 import Sage
from sage.misc.preparser import preparse_file

from sage.dsage.database.job import Job, expand_job
from sage.dsage.misc.hostinfo import ClassicHostInfo
from sage.dsage.errors.exceptions import NoJobException
from sage.dsage.twisted.pb import PBClientFactory
from sage.dsage.misc.constants import DELIMITER
from sage.dsage.misc.constants import DSAGE_DIR
from sage.dsage.misc.misc import random_str
from sage.dsage.misc.config import get_bool

START_MARKER = '___BEGIN___'
END_MARKER = '___END___'
LOG_PREFIX = "[Worker %s] "

def unpickle(pickled_job):
    return cPickle.loads(zlib.decompress(pickled_job))

class Worker(object):
    """
    This class represents a worker object that does the actual calculation.

    Parameters:
    remoteobj -- reference to the remote PB server

    """

    def __init__(self, remoteobj, id, log_level=0, delay=5.0):
        self.remoteobj = remoteobj
        self.id = id
        self.free = True
        self.job = None
        self.log_level = log_level
        self.delay = delay
        self.checker_task = task.LoopingCall(self.check_work)
        self.checker_timeout = 0.5
        self.got_output = False
        self.job_start_time = None
        self.sleep_time = 5.0
        self.start()

        # import some basic modules into our Sage() instance
        self.sage.eval('import time')
        self.sage.eval('import os')

    def _catch_failure(self, failure):
        log.msg("Error: ", failure.getErrorMessage())
        log.msg("Traceback: ", failure.printTraceback())

    def get_job(self):
        try:
            if self.log_level > 3:
                log.msg(LOG_PREFIX % self.id +  'Getting job...')
            d = self.remoteobj.callRemote('get_job')
        except Exception, msg:
            log.msg(msg)
            log.msg(LOG_PREFIX % self.id +  'Disconnected...')
            reactor.callLater(self.delay, self.get_job)
            return
        d.addCallback(self.gotJob)
        d.addErrback(self.noJob)

        return d

    def gotJob(self, jdict):
        """
        gotJob is a callback for the remoteobj's get_job method.

        Parameters:
        job -- Job object returned by remote's 'get_job' method

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
            self.sleep_time = 5.0 # Sets it to the default
            self.doJob(self.job)
        except Exception, msg:
            log.msg(msg)
            self.report_failure(msg)
            self.restart()

    def job_done(self, output, result, completed):
        """
        job_done is a callback for doJob.  Called when a job completes.

        Parameters:
        output -- the output of the command
        result -- the result of processing the job, a pickled object
        completed -- whether or not the job is completely finished (bool)

        """

        try:
            d = self.remoteobj.callRemote('job_done',
                                          self.job.job_id,
                                          output,
                                          result,
                                          completed)
        except Exception, msg:
            log.msg(msg)
            log.msg('[Worker: %s, job_done] Disconnected, reconnecting in %s'\
                    % (self.id, self.delay))
            reactor.callLater(self.delay, self.job_done,
                              output, result, completed)
            d = defer.Deferred()
            d.errback(error.ConnectionLost())

            return d

        if completed:
            self.restart()

        return d

    def noJob(self, failure):
        """
        noJob is an errback that catches the NoJobException.

        Parameters:
        failure -- a twisted.python.failure object (twisted.python.failure)

        """

        if failure.check(NoJobException):
            if self.sleep_time >= 60:
                self.sleep_time = 60
            else:
                self.sleep_time = self.sleep_time * 2
            if self.log_level > 3:
                msg = 'Sleeping for %s seconds' % self.sleep_time
                log.msg(LOG_PREFIX % self.id + msg)
            reactor.callLater(self.sleep_time, self.get_job)
        else:
            log.msg("Error: ", failure.getErrorMessage())
            log.msg("Traceback: ", failure.printTraceback())

    def setup_tmp_dir(self, job):
        """
        Creates the temporary directory for the worker.

        """

        cur_dir = os.getcwd() # keep a reference to the current directory
        tmp_dir = os.path.join(DSAGE_DIR, 'tmp_worker_files')
        tmp_job_dir = os.path.join(tmp_dir, job.job_id)
        if not os.path.isdir(tmp_dir):
            os.mkdir(tmp_dir)
        if not os.path.isdir(tmp_job_dir):
            os.mkdir(tmp_job_dir)
        os.chdir(tmp_job_dir)
        self.sage.eval("os.chdir('%s')" % tmp_job_dir)

        return tmp_job_dir

    def extract_job_data(self, job):
        """
        Extracts all the data that is in a job object.

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
                        f = open(var, 'wb')
                        f.write(data)
                        f.close()
                        if self.log_level > 2:
                            msg = 'Extracted %s' % f
                            log.msg(LOG_PREFIX % self.id + msg)
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

        """

        parsed_file = preparse_file(job.code, magic=False,
                                    do_time=False, ignore_prompts=False)

        job_filename = str(job.name) + '.py'
        job_file = open(job_filename, 'w')
        timeout = job.timeout
        BEGIN = "print '%s'\n\n" % (START_MARKER)
        END = "print '%s'\n\n" % (END_MARKER)
        SAVE_RESULT = """try:
    save(DSAGE_RESULT, 'result.sobj', compress=True)
except:
    save('No DSAGE_RESULT', 'result.sobj', compress=True)"""
        job_file.write("alarm(%s)\n\n" % (timeout))
        job_file.write(BEGIN)
        job_file.write(parsed_file)
        job_file.write("\n\n")
        job_file.write(END)
        job_file.write(SAVE_RESULT)
        job_file.close()
        if self.log_level > 2:
            log.msg('[Worker: %s] Wrote job file. ' % (self.id))

        return job_filename

    def doJob(self, job):
        """
        doJob is the method that drives the execution of a job.

        Parameters:
        job -- a Job object (dsage.database.Job)

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
        self.extract_job_data(job)

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
        try:
            os.chdir(self.tmp_job_dir)
            foo, output, new = self.sage._so_far(wait=0.1)
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
            self.free = True
            self.reset_checker()
        else:
            result = cPickle.dumps('Job not done yet.', 2)
        if self.check_failure(new):
            self.report_failure(new)
            self.restart()
            return
        sanitized_output = self.clean_output(new)
        if sanitized_output == '' and not done:
            self.increase_checker_task_timeout()
        else:
            d = self.job_done(sanitized_output, result, done)
            d.addErrback(self._catch_failure)

    def report_failure(self, failure):
        """
        Reports failure of a job.

        """

        import shutil
        shutil.move(self.tmp_job_dir, self.tmp_job_dir + '_failed')

        msg = 'Job %s failed!' % (self.job.job_id)
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

        self.checker_timeout = self.checker_timeout * 2
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
        try:
            pid = self.sage.pid()
            cmd = 'kill -9 -%s'%pid
            os.system(cmd)
            self.start()
        except Exception, msg:
            self.start()

    def stop(self, hard_reset=False):
        """
        Stops the current worker and resets it's internal state.

        """

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
                        if self.log_level > 3:
                            msg = 'Interrupting SAGE (try %s)' % i
                            log.msg(LOG_PREFIX % self.id + msg)
            except Exception, msg:
                success = False
                log.msg(msg)
                log.msg(LOG_PREFIX % self.id + "Performing hard reset.")

            if not success:
                self.kill_sage()
                self.start()
            else:
                self.sage.reset()
        self.free = True
        self.job = None

    def start(self):
        """
        Starts a new worker if it does not exist already.

        """

        if not hasattr(self, 'sage'):
            if self.log_level > 3:
                logfile = DSAGE_DIR + '/%s-pexpect.log' % self.id
                self.sage = Sage(maxread=1, logfile=logfile)
            else:
                self.sage = Sage(maxread=1)
            try:
                self.sage._start(block_during_init=True)
            except RuntimeError, msg: # Could not start SAGE
                print msg
                print 'Failed to start a worker, probably Expect issues.'
                reactor.stop()
                sys.exit(-1)

        self.get_job()

    def restart(self):
        """
        Restarts the current worker.

        """

        try:
            delta = datetime.datetime.now() - self.job_start_time
            if delta.seconds >= (60*5): # more than 5 minutes, do a hard reset
                self.stop(hard_reset=True)
            else:
                self.stop()
        except TypeError:
            self.stop(hard_reset=True)
        self.job_start_time = None
        self.start()
        self.reset_checker()
        log.msg('[Worker: %s] Restarting...' % (self.id))

class Monitor(object):
    """
    This class represents a monitor that controls workers.

    It monitors the workers and checks on their status

    Parameters:
    hostname -- the hostname of the server we want to connect to (str)
    port -- the port of the server we want to connect to (int)

    """

    def __init__(self, server='localhost', port=8081,
                 username=getuser(),
                 ssl=True,
                 workers=2,
                 anonymous=False,
                 priority=20,
                 delay=5.0,
                 log_level=0,
                 log_file=os.path.join(DSAGE_DIR, 'worker.log'),
                 pubkey_file=None,
                 privkey_file=None):
        self.server = server
        self.port = port
        self.username = username
        self.ssl = ssl
        self.workers = workers
        self.anonymous = anonymous
        self.priority = priority
        self.delay = delay
        self.log_level = log_level
        self.log_file = log_file
        self.pubkey_file = pubkey_file
        self.privkey_file = privkey_file

        self.remoteobj = None
        self.connected = False
        self.reconnecting = False
        self.worker_pool = None
        self.sleep_time = 1.0

        self.host_info = ClassicHostInfo().host_info
        try:
            uuid_file = os.path.join(DSAGE_DIR, 'worker.uuid')
            uuid_ = open(uuid_file).readlines()[0]
        except Exception, e:
            uuid_ = str(uuid.uuid1())
            f = open(uuid_file, 'w')
            f.write(uuid_)
            f.close()
        self.host_info['uuid'] = uuid_
        self.host_info['workers'] = self.workers

        self._startLogging(self.log_file)

        try:
            os.nice(self.priority)
        except OSError, msg:
            log.msg('Error setting priority: %s' % (self.priority))
            pass
        if not self.anonymous:
            from twisted.cred import credentials
            from twisted.conch.ssh import keys
            self.DATA =  random_str(500)
            # public key authentication information
            self.pubkey_str =keys.getPublicKeyString(self.pubkey_file)
            # try getting the private key object without a passphrase first
            try:
                self.priv_key = keys.getPrivateKeyObject(self.privkey_file)
            except keys.BadKeyError:
                pphrase = self._getpassphrase()
                self.priv_key = keys.getPrivateKeyObject(self.privkey_file,
                                                         pphrase)
            self.pub_key = keys.getPublicKeyObject(self.pubkey_str)
            self.alg_name = 'rsa'
            self.blob = keys.makePublicKeyBlob(self.pub_key)
            self.data = self.DATA
            self.signature = keys.signData(self.priv_key, self.data)
            self.creds = credentials.SSHPrivateKey(self.username,
                                                   self.alg_name,
                                                   self.blob,
                                                   self.data,
                                                   self.signature)

    def _startLogging(self, log_file):
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
        self.remoteobj = remoteobj
        self.remoteobj.notifyOnDisconnect(self._disconnected)
        self.connected = True
        self.reconnecting = False

        if self.worker_pool == None: # Only pool workers the first time
            self.pool_workers(self.remoteobj)
        else:
            for worker in self.worker_pool:
                worker.remoteobj = self.remoteobj # Update workers

    def _disconnected(self, remoteobj):
        log.msg('Lost connection to the server.')
        self.connected = False
        self._retryConnect()

    def _got_killed_jobs(self, killed_jobs):
        if killed_jobs == None:
            return
        # reconstruct the Job objects from the jdicts
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
        log.msg('[Monitor] Disconnected, reconnecting in %s' % self.delay)
        if not self.connected:
            reactor.callLater(self.delay, self.connect)

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

        factory = pb.PBClientFactory()
        self.factory = PBClientFactory()
        try:
            if self.ssl:
                from twisted.internet import ssl
                contextFactory = ssl.ClientContextFactory()
                reactor.connectSSL(self.server, self.port,
                                   self.factory, contextFactory)
            else:
                reactor.connectTCP(self.server, self.port, self.factory)
        except Exception, msg:
            self._retryConnect()

        log.msg(DELIMITER)
        log.msg('DSAGE Worker')
        log.msg('Started with PID: %s' % (os.getpid()))
        log.msg('Connecting to %s:%s' % (self.server, self.port))
        if self.ssl:
            log.msg('Using SSL: True')
        else:
            log.msg('Using SSL: False')
        log.msg(DELIMITER)

        if not self.anonymous:
            log.msg('Connecting as authenticated worker...\n')
            d = self.factory.login(self.creds,
                                   (pb.Referenceable(), self.host_info))
        else:
            log.msg('Connecting as anonymous worker...\n')
            d = self.factory.login('Anonymous',
                                   (pb.Referenceable(), self.host_info))
        d.addCallback(self._connected)
        d.addErrback(self._catchConnectionFailure)

        return d

    def pool_workers(self, remoteobj):
        """
        pool_workers creates as many workers as specified in worker.conf.

        """

        log.msg('[Monitor] Starting %s workers...' % (self.workers))
        self.worker_pool = [Worker(remoteobj, x, self.log_level, self.delay)
                            for x in range(self.workers)]

    def check_killed_jobs(self):
        """
        check_killed_jobs retrieves a list of killed job ids.

        """

        if not self.connected:
            return

        killed_jobs = self.remoteobj.callRemote('get_killed_jobs_list')
        killed_jobs.addCallback(self._got_killed_jobs)
        killed_jobs.addErrback(self._catch_failure)

    def start_looping_calls(self):
        """
        start_looping_calls prepares and starts our periodic checking methods.

        """
        #
        # self.check_output_timeout = 1
        # self.tsk1 = task.LoopingCall(self.check_output)
        # self.tsk1.start(self.check_output_timeout, now=False)

        interval = 5.0
        self.tsk2 = task.LoopingCall(self.check_killed_jobs)
        self.tsk2.start(interval, now=False)

    def stop_looping_calls(self):
        """
        stops the looping calls.

        """

        # self.tsk1.stop()
        self.tsk2.stop()

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
    parser.add_option('-d', '--delay',
                      dest='delay',
                      type='float',
                      default=5,
                      help='delay before checking for new job. default=5')
    parser.add_option('-a', '--anonymous',
                      dest='anonymous',
                      default=False,
                      action='store_true',
                      help='Connect as anonymous worker. default=False')
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
                      default=True,
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
    SSL = get_bool(options.ssl)
    monitor = Monitor(server=options.server,
                      port=options.port,
                      username=options.username,
                      ssl=SSL,
                      workers=options.workers,
                      anonymous=options.anonymous,
                      priority=options.priority,
                      delay=options.delay,
                      log_file=options.logfile,
                      log_level=options.loglevel,
                      pubkey_file=options.pubkey_file,
                      privkey_file=options.privkey_file)
    monitor.connect()
    monitor.start_looping_calls()
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