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
import ConfigParser
import uuid
import cPickle
import zlib

from twisted.spread import pb
from twisted.internet import reactor, defer, error, task
from twisted.python import log

from sage.interfaces.sage0 import Sage
from sage.misc.preparser import preparse_file

from sage.dsage.database.job import Job, expand_job
from sage.dsage.misc.hostinfo import HostInfo, ClassicHostInfo
from sage.dsage.errors.exceptions import NoJobException
from sage.dsage.twisted.pb import PBClientFactory
from sage.dsage.misc.constants import delimiter as DELIMITER

pb.setUnjellyableForClass(HostInfo, HostInfo)

DSAGE_DIR = os.path.join(os.getenv('DOT_SAGE'), 'dsage')

# Begin reading configuration
try:
    CONF_FILE = os.path.join(DSAGE_DIR, 'worker.conf')
    CONFIG = ConfigParser.ConfigParser()
    CONFIG.read(CONF_FILE)

    LOG_FILE = CONFIG.get('log', 'log_file')
    LOG_LEVEL = CONFIG.getint('log','log_level')
    SSL = CONFIG.getint('ssl', 'ssl')
    WORKERS = CONFIG.getint('general', 'workers')
    SERVER = CONFIG.get('general', 'server')
    PORT = CONFIG.getint('general', 'port')
    DELAY = CONFIG.getint('general', 'delay')
    NICE_LEVEL = CONFIG.getint('general', 'nice_level')
    ANONYMOUS = CONFIG.getboolean('general', 'anonymous')
except Exception, msg:
    print msg
    print "Error reading %s, please fix manually or run dsage.setup()" % \
    CONF_FILE
    sys.exit(-1)
# End reading configuration

# OUTPUT MARKERS shared by Worker and Monitor
START_MARKER = '___BEGIN___'
END_MARKER = '___END___'

def unpickle(pickled_job):
    return cPickle.loads(zlib.decompress(pickled_job))

class Worker(object):
    r"""
    This class represents a worker object that does the actual calculation.

    Parameters:
    remoteobj -- reference to the remote PB server

    """

    def __init__(self, remoteobj, id):
        self.remoteobj = remoteobj
        self.id = id
        self.free = True
        self.job = None

        if LOG_LEVEL > 3:
            self.sage = Sage(logfile=DSAGE_DIR + '/%s-pexpect.log'\
                             % self.id)
        else:
            self.sage = Sage()

        # import some basic modules into our Sage() instance
        self.sage.eval('import time')
        self.sage.eval('import sys')
        self.sage.eval('import os')

        # Initialize getting of jobs
        self.get_job()

    def get_job(self):
        try:
            if LOG_LEVEL > 3:
                log.msg('Worker %s: Getting job...' % (self.id))
            d = self.remoteobj.callRemote('get_job')
        except Exception, msg:
            log.msg(msg)
            log.msg('[Worker: %s, get_job] Disconnected from remote server.'\
                    % self.id)
            reactor.callLater(DELAY, self.get_job)
            return
        d.addCallback(self.gotJob)
        d.addErrback(self.noJob)

        return d

    def gotJob(self, jdict):
        r"""
        gotJob is a callback for the remoteobj's get_job method.

        Parameters:
        job -- Job object returned by remote's 'get_job' method

        """

        if LOG_LEVEL > 3:
            log.msg('[Worker %s, gotJob] %s' % (self.id, jdict))

        self.job = expand_job(jdict)

        if not isinstance(self.job, Job):
            raise NoJobException

        log.msg('[Worker: %s] Got job (%s, %s)' % (self.id,
                                                   self.job.name,
                                                   self.job.id))
        try:
            self.doJob(self.job)
        except Exception, e:
            print e
            raise

    def job_done(self, output, result, completed):
        r"""
        job_done is a callback for doJob.  Called when a job completes.

        Parameters:
        output -- the output of the command
        result -- the result of processing the job, a pickled object
        completed -- whether or not the job is completely finished (bool)

        """

        try:
            d = self.remoteobj.callRemote('job_done',
                                          self.job.id,
                                          output,
                                          result,
                                          completed)
        except Exception, msg:
            log.msg(msg)
            log.msg('[Worker: %s, job_done] Disconnected, reconnecting in %s'\
                    % (self.id, DELAY))
            reactor.callLater(DELAY, self.job_done, output, result, completed)
            d = defer.Deferred()
            d.errback(error.ConnectionLost())
            return d

        if completed:
            self.restart()

        return d

    def noJob(self, failure):
        # TODO: Probably do not need this errback, look into consolidating
        # with failedJob
        r"""
        noJob is an errback that catches the NoJobException.

        Parameters:
        failure -- a twisted.python.failure object (twisted.python.failure)

        """

        sleep_time = 5.0
        if failure.check(NoJobException):
            if LOG_LEVEL > 3:
                log.msg('[Worker %s, noJob] Sleeping for %s seconds\
                ' % (self.id, sleep_time))
            reactor.callLater(5.0, self.get_job)
        else:
            print "Error: ", failure.getErrorMessage()
            print "Traceback: ", failure.printTraceback()

    def setup_tmp_dir(self, job):
        # change to a temporary directory
        cur_dir = os.getcwd() # keep a reference to the current directory
        tmp_dir = os.path.join(DSAGE_DIR, 'tmp_worker_files')
        tmp_job_dir = os.path.join(tmp_dir, job.id)
        if not os.path.isdir(tmp_dir):
            os.mkdir(tmp_dir)
        os.mkdir(tmp_job_dir)
        os.chdir(tmp_job_dir)
        self.sage.eval("os.chdir('%s')" % tmp_job_dir)

        return tmp_job_dir

    def extract_job_data(self, job):
        r"""
        Extracts all the data that is in a job object.

        """

        if isinstance(job.data, list):
            if LOG_LEVEL > 2:
                log.msg('Extracting job data...')
            for var, data, kind in job.data:
                # Uncompress data
                try:
                    data = zlib.decompress(data)
                except Exception, msg:
                    log.msg(msg)
                    continue
                if kind == 'file':
                    # Write out files to current dir
                    f = open(var, 'wb')
                    f.write(data)
                    if LOG_LEVEL > 2:
                        log.msg('[Worker: %s] Extracted %s. ' % (self.id, f))
                if kind == 'object':
                    # Load object into the SAGE worker
                    fname = var + '.sobj'
                    if LOG_LEVEL > 3:
                        log.msg('Object to be loaded: %s' % fname)
                    f = open(fname, 'wb')
                    f.write(data)
                    f.close()
                    self.sage.eval("%s = load('%s')" % (var, fname))
                    if LOG_LEVEL > 2:
                        log.msg('[Worker: %s] Loaded %s' % (self.id, fname))

    def write_job_file(self, job):
        r"""
        Writes out the job file to be executed to disk.

        """
        parsed_file = preparse_file(job.code, magic=False,
                                    do_time=False, ignore_prompts=False)

        job_filename = str(job.name) + '.py'
        job_file = open(job_filename, 'w')
        BEGIN = "print '%s'\n\n" % (START_MARKER)
        END = "print '%s'\n\n" % (END_MARKER)
        job_file.write(BEGIN)
        job_file.write(parsed_file)
        job_file.write("\n\n")
        job_file.write(END)
        job_file.close()
        if LOG_LEVEL > 2:
            log.msg('[Worker: %s] Wrote job file. ' % (self.id))

        return job_filename

    def doJob(self, job):
        r"""
        doJob is the method that drives the execution of a job.

        Parameters:
        job -- a Job object (dsage.database.Job)

        """

        if LOG_LEVEL > 3:
            log.msg('[Worker %s, doJob] Beginning job execution...' % (self.id))

        self.free = False
        d = defer.Deferred()

        self.tmp_job_dir = self.setup_tmp_dir(job)
        self.extract_job_data(job)

        job_filename = self.write_job_file(job)

        f = os.path.join(self.tmp_job_dir, job_filename)
        self.sage._send("execfile('%s')" % (f))
        if LOG_LEVEL > 2:
            log.msg('[Worker: %s] File to execute: %s' % (self.id, f))

        d.callback(True)

        return d

    def stop(self):
        r"""
        Stops the current worker and resets it's internal state.

        """

        # This quits the current running calculation
        self.sage._expect.sendline(chr(3))  # send ctrl-c
        self.sage._expect.expect(self.sage._prompt)
        self.sage._expect.expect(self.sage._prompt)

        self.sage.reset()
        self.free = True
        self.job = None

    def start(self):
        r"""
        Starts a new worker if it does not exist already.

        """
        if self.sage is None:
            if LOG_LEVEL > 3:
                self.sage = Sage(logfile=DSAGE_DIR + '/%s-pexpect.out' % self.id)
            else:
                self.sage = Sage()
        self.get_job()

    def restart(self):
        r"""
        Restarts the current worker.

        """

        log.msg('[Worker: %s] Restarting...' % (self.id))
        self.stop()
        self.start()

class Monitor(object):
    r"""
    This class represents a monitor that controls workers.

    It monitors the workers and checks on their status

    Parameters:
    hostname -- the hostname of the server we want to connect to (str)
    port -- the port of the server we want to connect to (int)

    """

    def __init__(self, hostname, port):
        if hostname is None:
            hostname = SERVER
        self.hostname = hostname
        if port is None:
            port = PORT
        self.port = port
        self.remoteobj = None
        self.connected = False
        self.reconnecting = False
        self.workers = None

        # Start twisted logging facility
        self._startLogging(LOG_FILE)

        if len(CONFIG.get('uuid', 'id')) != 36:
            CONFIG.set('uuid', 'id', str(uuid.uuid1()))
            f = open(CONF_FILE, 'w')
            CONFIG.write(f)
        self.uuid = CONFIG.get('uuid', 'id')

        self.host_info = ClassicHostInfo().host_info
        self.host_info['uuid'] = self.uuid
        self.host_info['workers'] = WORKERS

        if not ANONYMOUS:
            from twisted.cred import credentials
            from twisted.conch.ssh import keys
            self._get_auth_info()
            # public key authentication information
            self.pubkey_str =keys.getPublicKeyString(filename=self.pubkey_file)
            # try getting the private key object without a passphrase first
            try:
                self.priv_key = keys.getPrivateKeyObject(filename=self.privkey_file)
            except keys.BadKeyError:
                pphrase = self._getpassphrase()
                self.priv_key = keys.getPrivateKeyObject(filename=self.privkey_file, passphrase=pphrase)
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
        else:
            print "Logging to file: ", log_file
            server_log = open(log_file, 'a')
            log.startLogging(server_log)

    def _get_auth_info(self):
        import random
        self.DATA =  ''.join([chr(i) for i in [random.randint(65, 123) for n in
                        range(500)]])
        self.DSAGE_DIR = os.path.join(os.getenv('DOT_SAGE'), 'dsage')
        # Begin reading configuration
        try:
            conf_file = os.path.join(self.DSAGE_DIR, 'client.conf')
            config = ConfigParser.ConfigParser()
            config.read(conf_file)

            self.port = config.getint('general', 'port')
            self.username = config.get('auth', 'username')
            self.privkey_file = os.path.expanduser(config.get('auth', 'privkey_file'))
            self.pubkey_file = os.path.expanduser(config.get('auth', 'pubkey_file'))
        except Exception, msg:
            print msg
            raise

    def _getpassphrase(self):
        import getpass
        passphrase = getpass.getpass('Passphrase (Hit enter for None): ')

        return passphrase

    def _connected(self, remoteobj):
        self.remoteobj = remoteobj
        self.remoteobj.notifyOnDisconnect(self._disconnected)
        self.connected = True
        self.reconnecting = False

        if self.workers == None: # Only pool workers the first time
            self.pool_workers(self.remoteobj)
        else:
            for worker in self.workers:
                worker.remoteobj = self.remoteobj # Update workers
        # self.submit_host_info()

    def _disconnected(self, remoteobj):
        log.msg('Lost connection to the server.')
        self.connected = False
        self._retryConnect()

    def _got_killed_jobs(self, killed_jobs):
        if killed_jobs == None:
            return
        # reconstruct the Job objects from the jdicts
        killed_jobs = [expand_job(jdict) for jdict in killed_jobs]
        for worker in self.workers:
            if worker.job is None:
                continue
            if worker.free:
                continue
            for job in killed_jobs:
                if job is None or worker.job is None:
                    continue
                if worker.job.id == job.id:
                    log.msg('[Worker: %s] Processing a killed job, \
                            restarting...' % worker.id)
                    worker.restart()

    def _retryConnect(self):
        log.msg('[Monitor] Disconnected, reconnecting in %s' % DELAY)
        reactor.callLater(DELAY, self.connect)

    def _catchConnectionFailure(self, failure):
        # If we lost the connection to the server somehow
        # if failure.check(error.ConnectionRefusedError,
        #                 error.ConnectionLost,
        #                 pb.DeadReferenceError):

        self.connected = False
        self._retryConnect()

        log.msg("Error: ", failure.getErrorMessage())
        log.msg("Traceback: ", failure.printTraceback())

    def _catch_failure(self, failure):
        log.msg("Error: ", failure.getErrorMessage())
        log.msg("Traceback: ", failure.printTraceback())

    def connect(self):
        r"""
        This method connects the monitor to a remote PB server.

        """
        if self.connected: # Don't connect multiple times
            return

        factory = pb.PBClientFactory()

        log.msg(DELIMITER)
        log.msg('DSAGE Worker')
        log.msg('Connecting to %s:%s' % (self.hostname, self.port))
        log.msg(DELIMITER)

        self.factory = PBClientFactory()
        if SSL == 1:
            from twisted.internet import ssl
            contextFactory = ssl.ClientContextFactory()
            reactor.connectSSL(self.hostname, self.port,
                               self.factory, contextFactory)
        else:
            reactor.connectTCP(self.hostname, self.port, self.factory)

        if not ANONYMOUS:
            log.msg('Connecting as authenticated worker...\n')
            d = self.factory.login(self.creds, (pb.Referenceable(), self.host_info))
        else:
            log.msg('Connecting as anonymous worker...\n')
            d = self.factory.login('Anonymous', (pb.Referenceable(), self.host_info))
        d.addCallback(self._connected)
        d.addErrback(self._catchConnectionFailure)

        return d

    def pool_workers(self, remoteobj):
        r"""
        pool_workers creates as many workers as specified in worker.conf.

        """

        self.workers = [Worker(remoteobj, x) for x in range(WORKERS)]
        log.msg('[Monitor] Initialized ' + str(WORKERS) + ' workers.')

    def check_output(self):
        r"""
        check_output periodically polls workers for new output.

        This figures out whether or not there is anything new output that we
        should submit to the server.

        """

        if self.workers == None:
            return

        for worker in self.workers:
            if worker.job == None:
                continue
            if worker.free == True:
                continue

            if LOG_LEVEL > 1:
                log.msg('[Monitor] Checking for job output')
            try:
                done, output, new = worker.sage._so_far()
            except Exception, msg:
                log.msg('Exception raised when checking output.')
                log.msg(msg)
                continue
            if new == '' or new.isspace():
                continue
            if done:
                worker.free = True
                # Checks to see if the job created a result var
                # Do this multiple times because of Expect weirdness
                sobj = worker.sage.get('DSAGE_RESULT')
                if sobj == '' or sobj.isspace():
                    sobj = worker.sage.get('DSAGE_RESULT')
                    if sobj == '' or sobj.isspace():
                        sobj = worker.sage.get('DSAGE_RESULT')
                    else:
                        if LOG_LEVEL > 1:
                            log.msg('Got DSAGE_RESULT second time')

                if 'Error: name \'DSAGE_RESULT\' is not defined' in sobj:
                    if LOG_LEVEL > 1:
                        log.msg('DSAGE_RESULT does not exist')
                    result = cPickle.dumps('No result saved.', 2)
                else:
                    worker.sage.eval("save(DSAGE_RESULT, 'result')")
                    os.chdir(worker.tmp_job_dir)
                    try:
                        result = open('result.sobj', 'rb').read()
                    except Exception, msg:
                        if LOG_LEVEL > 1:
                            log.msg(msg)
                        result = cPickle.dumps('Error in reading result.', 2)
                log.msg("Job '%s' finished" % worker.job.id)
            else:
                result = cPickle.dumps('Job not done yet.', 2)

            sanitized_output = self.clean_output(new)

            if self.check_failure(sanitized_output):
                s = ['[Monitor] Error in result for ',
                     'job %s %s done by ' % (worker.job.name, worker.job.id),
                     'Worker %s' % (worker.id)
                     ]
                log.msg(''.join(s))
                log.msg('[Monitor] Traceback: \n%s' % sanitized_output)
                d = self.remoteobj.callRemote('job_failed', worker.job.id, sanitized_output)
                d.addErrback(self._catch_failure)
                continue

            d = worker.job_done(sanitized_output, result, done)
            d.addErrback(self._catchConnectionFailure)

    def check_failure(self, sage_output):
        r"""
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

    def check_killed_jobs(self):
        r"""
        check_killed_jobs retrieves a list of killed job ids.

        """

        if not self.connected:
            return

        killed_jobs = self.remoteobj.callRemote('get_killed_jobs_list')
        killed_jobs.addCallback(self._got_killed_jobs)
        killed_jobs.addErrback(self._catch_failure)

    def clean_output(self, sage_output):
        r"""
        clean_output attempts to clean up the output string from sage.

        """

        # log.msg("Before cleaning output: ", sage_output)
        begin = sage_output.find(START_MARKER)
        if begin != -1:
            begin += len(START_MARKER)
        else:
            begin = 0
        end = sage_output.find(END_MARKER)
        if end != -1:
            end -= 1
        else:
            end = len(sage_output)
        output = sage_output[begin:end]
        output = output.strip()
        output = output.replace('\r', '')

        # log.msg("After cleaning output: ", output)
        return output

    def start_looping_calls(self):
        r"""
        stop_looping_calls prepares and starts our periodic checking methods.

        """

        # submits the output to the server
        self.tsk1 = task.LoopingCall(self.check_output)
        self.tsk1.start(0.1, now=False)
        # checks for killed jobs
        self.tsk2 = task.LoopingCall(self.check_killed_jobs)
        self.tsk2.start(5.0, now=False)

    def stop_looping_calls(self):
        r"""
        Stops the looping calls.

        """
        self.tsk.stop()
        self.tsk1.stop()
        self.tsk2.stop()

def main():
    r"""
    argv[1] == hostname
    argv[2] == port

    """

    if len(sys.argv) == 3:
        hostname, port = sys.argv[1:3]
        try:
            port = int(port)
        except:
            port = None
        if hostname == 'None':
            hostname = None
        else:
            try:
                hostname = str(hostname)
            except Exception, msg:
                print msg
                hostname = None
    else:
        hostname = port = None

    monitor = Monitor(hostname, port)

    monitor.connect()
    monitor.start_looping_calls()

    try:
        reactor.run()
    except:
        sys.exit()

if __name__ == '__main__':
    main()

