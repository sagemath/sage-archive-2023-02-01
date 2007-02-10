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
############################################################################

import sys
import os
import random
import glob
import ConfigParser
import copy
import cPickle
import zlib
import time
import threading

from twisted.python import threadable
from twisted.spread import pb
from twisted.internet import reactor, defer, error, task
from twisted.cred import credentials
from twisted.conch.ssh import keys

from sage.dsage.database.job import Job
from sage.dsage.pb.dsage_pb import ClientPBClientFactory
from sage.dsage.errors.exceptions import NoJobException, NotConnectedException

# This is a randomly generated string we will use to as the signature to
# sign
DATA =  ''.join([chr(i) for i in [random.randint(65, 123) for n in
                range(500)]])

DSAGE_DIR = os.path.join(os.getenv('DOT_SAGE'), 'dsage')
# Begin reading configuration
try:
    conf_file = os.path.join(DSAGE_DIR, 'client.conf')
    config = ConfigParser.ConfigParser()
    config.read(conf_file)

    LOG_FILE = config.get('log', 'log_file')
    LOG_LEVEL = config.getint('log', 'log_level')
    SSL = config.getint('ssl', 'ssl')
    USERNAME = config.get('auth', 'username')
    PRIVKEY_FILE = os.path.expanduser(config.get('auth', 'privkey_file'))
    PUBKEY_FILE = os.path.expanduser(config.get('auth', 'pubkey_file'))
    SERVER = config.get('general', 'server')
    PORT = config.getint('general', 'port')
except:
    print 'Error reading %s, please fix manually or run dsage.setup()' % \
    (conf_file)
    sys.exit(-1)
# End reading configuration

class DSageThread(threading.Thread):
    def run(self):
        reactor.run(installSignalHandlers=False)

class DSage(object):
    r"""
    This object represents a connection to the server

       Parameters:
       hostname -- hostname of the DSAGE server (str)
       port -- port of the server (int)
       username -- username stored on the server (str)
       pubkey_file -- file that stores the users public key
       privkey_file -- file that stores the users private key

    """
    def __init__(self, server=None, port=8081, username=None,
                 pubkey_file=None, privkey_file=None):

        # We will read the values in from the conf file first and let the
        # user override the values stored in the conf file by keyword
        # parameters

        if server is None:
            self.server = SERVER
        else:
            self.server = server

        self.port = PORT
        self.username = USERNAME
        self.pubkey_file = PUBKEY_FILE
        self.privkey_file = PRIVKEY_FILE

        self.remoteobj = None
        self.result = None

        # public key authentication information
        self.public_key_string = keys.getPublicKeyString(
                                                    filename=self.pubkey_file)
        self.private_key = keys.getPrivateKeyObject(filename=self.privkey_file)
        self.public_key = keys.getPublicKeyObject(self.public_key_string)
        self.alg_name = 'rsa'
        self.blob = keys.makePublicKeyBlob(self.public_key)
        self.data = DATA
        self.signature = keys.signData(self.private_key, self.data)
        self.creds = credentials.SSHPrivateKey(self.username,
                                               self.alg_name,
                                               self.blob,
                                               self.data,
                                               self.signature)

        self.jobs = []

        self.connect()

    def __str__(self):
        self.check_connected()
        self.info_str = 'Connected to: ' \
                    + self.server + ':' + str(self.port)
        return self.info_str + '\r'

    def __call__(self, cmd, job_name=None):
        cmd = ['ans = %s\n' % (cmd),
               'print ans\n'
               "save(ans, 'ans')\n"
               "DSAGE_RESULT = 'ans.sobj'\n"]

        return self.eval(''.join(cmd), job_name)

    def __getstate__(self):
        d = copy.copy(self.__dict__)
        d['remoteobj'] = None
        return d

    def _catchFailure(self, failure):
        if failure.check(error.ConnectionRefusedError):
            print 'Remote server refused our connection. \r'
            return
        print "Error: ", failure.getErrorMessage()
        print "Traceback: ", failure.printTraceback()

    def _connected(self, remoteobj):
        if LOG_LEVEL > 0:
            print 'Connected to remote server.\r'
        self.remoteobj = remoteobj
        self.remoteobj.notifyOnDisconnect(self._disconnected)

    def _disconnected(self, remoteobj):
        print 'Lost connection to the server.'
        self.info_str = 'Not connected.'

    def _gotMyJobs(self, jobs, job_name):
        if jobs == None:
            raise NoJobException
        if job_name:
            return [JobWrapper(self.remoteobj, job)
                    for job in jobs if job.name == job_name]

    def _killedJob(self, jobID):
        if jobID:
            if LOG_LEVEL > 2:
                print str(jobID) + ' was successfully killed.'

    def restore(self, remoteobj):
        r"""
        This method restores a connection to the server.
        """
        self.remoteobj = remoteobj

    def connect(self):
        r"""
        This methods establishes the conection to the remote server.

        """

        # TODO: Send a useful 'mind' object with the login request!
        # factory = pb.PBClientFactory()
        factory = ClientPBClientFactory()

        if SSL == 1:
            from twisted.internet import ssl
            contextFactory = ssl.ClientContextFactory()
            reactor.connectSSL(self.server,
                               self.port,
                               factory,
                               contextFactory)
        else:
            reactor.connectTCP(self.server, self.port, factory)

        return factory.login(self.creds, None).addCallback(
                            self._connected).addErrback(
                            self._catchFailure)

    def disconnect(self):
        print 'Disconnecting from server.'
        self.remoteobj = None

    def eval(self, cmd, globals=None, job_name=None):
        r"""
        eval evaluates a command

        Parameters:
        cmd -- the sage command to be evaluated (str)
        globals -- a dict (see help for python's eval method)
        job_name -- an alphanumeric job name

        """

        self.check_connected()
        if not job_name or not isinstance(job_name, str):
            job_name = 'default_job'

        type = 'sage'

        job = Job(id=None, file=cmd, name=job_name,
                  author=self.username, type=type)

        wrapped_job = JobWrapper(self.remoteobj, job)
        if globals is not None:
            for k, v in globals.iteritems():
                job.attach(k, v)

        return wrapped_job

    def eval_file(self, fname, job_name):
        r"""
        eval_file allows you to evaluate the contents of an entire file.

        Parameters:
            fname -- file name of the file you wish to evaluate

        """

        self.check_connected()

        type = 'file'
        cmd = open(fname).read()
        # d = self.remoteobj.callRemote('getNextJobID')
        # d.addCallback(self._gotID, cmd, job_name, type)
        job = Job(id=None, file=cmd, name=job_name,
                  author=self.username, type=type)

        wrapped_job = JobWrapper(self.remoteobj, job)
        return wrapped_job

    def send_job(self, job):
        r"""
        Sends a Job object to the server.

        """

        if not isinstance(job, Job):
            raise TypeError
        wrapped_job = JobWrapper(self.remoteobj, job)
        # self.jobs.append(wrapped_job)
        return wrapped_job

    def _gotJobID(self, id, job):
        job.id = id
        job.author = self.username

        self.jobs.append(job)

        pickled_job = job.pickle()
        d = self.remoteobj.callRemote('submitJob', pickled_job)
        d.addErrback(self._catchFailure)
        # d.addCallback(self._submitted, job)

        return JobWrapper(self.remoteobj, job)

    def eval_dir(self, dir, job_name):
        self.check_connected()
        os.chdir(dir)
        files = glob.glob('*.spyx')
        deferreds = []
        for file in files:
            sage_cmd = open(file).readlines()
            d = self.remoteobj.callRemote('getNextJobID')
            d.addCallback(self._gotID, sage_cmd, job_name, file=True,
                          type='spyx')
            d.addErrback(self._catchFailure)
            deferreds.append(d)
        d_list = defer.DeferredList(deferreds)
        return d_list

    def kill(self, jobID):
        r"""
        Kills a job given the job id.

        Parameters:
        jobID -- job id

        """

        d = self.remoteobj.callRemote('killJob', jobID)
        d.addCallback(self._killedJob)
        d.addErrback(self._catchFailure)

    def get_my_jobs(self, is_active=False, job_name=None):
        r"""
        This method returns a list of jobs that belong to you.

        Parameters:
        is_active -- set to true to get only active jobs (bool)

        Use this method if you get disconnected from the server and wish to
        retrieve your old jobs back.

        """

        self.check_connected()

        d = self.remoteobj.callRemote('getJobsByAuthor',
                                      self.username,
                                      is_active,
                                      job_name)
        d.addCallback(self._gotMyJobs, job_name)
        d.addErrback(self._catchFailure)

        return d

    def cluster_speed(self):
        r"""
        Returns the speed of the cluster.

        """

        self.check_connected()

        return self.remoteobj.callRemote('getClusterSpeed')

    def check_connected(self):
        if self.remoteobj == None:
            raise NotConnectedException

    def checkResult(self, jobs):
        """Takes a list of jobs and checks the result."""


class JobWrapper(object):
    r"""
    Represents a remote job.

    Parameters:
        remoteobj -- the PB server's remoteobj
        job -- a Job object (job)

    """

    def __init__(self, remoteobj, job):
        self.remoteobj = remoteobj
        self._job = job

        # TODO Make this more complete
        self._update_job(job)
        self.worker_info = self._job.worker_info

        d = self.remoteobj.callRemote('getNextJobID')
        d.addCallback(self._gotID)
        d.addErrback(self._catchFailure)

        # hack to try and fetch a result after submitting the job
        self.syncJob_task = task.LoopingCall(self.syncJob)
        self.syncJob_task.start(2.0, now=True)
        # reactor.callLater(2.0, self.syncJob)
        # reactor.callLater(5.0, self.syncJob)
        # reactor.callLater(10.0, self.syncJob)

    def __repr__(self):
        if self._job.status == 'completed' and not self._job.output:
            return 'No output.'
        elif not self._job.output:
            return 'No output yet.'

        return self._job.output

    def __getstate__(self):
        d = copy.copy(self.__dict__)
        d['remoteobj'] = None
        d['syncJob_task'] = None
        return d

    def _update_job(self, job):
        for k, v in Job.__dict__.iteritems():
            if isinstance(v, property):
                setattr(self, k, getattr(job, k))
                #setattr(self, k ,self._job.jdic[k])

    def unpickle(self, pickled_job):
        return cPickle.loads(zlib.decompress(pickled_job))

    def wait(self):
        timeout = 0.1
        while self._job.result is None:
            reactor.iterate(timeout)

    def save(self, filename=None):
        if filename is None:
            filename = str(self._job.name)
        filename += '.sobj'
        f = open(filename, 'w')
        cPickle.dump(self, f, 2)
        return filename

    def _catchFailure(self, failure):
        if failure.check(pb.DeadReferenceError, error.ConnectionLost):
            print 'Disconnected from server.'
        else:
            print "Error: ", failure.getErrorMessage()
            print "Traceback: ", failure.printTraceback()

    def _gotJob(self, job):
        if job == None:
            return
        self._job = self.unpickle(job)
        self._update_job(self._job)

    def _gotID(self, id):
        self._job.id = id
        pickled_job = self._job.pickle()
        d = self.remoteobj.callRemote('submitJob', pickled_job)
        d.addErrback(self._catchFailure)

    def getJob(self):
        if self.remoteobj == None:
            raise NotConnectedException
        d = self.remoteobj.callRemote('getJobByID', self._job.id)
        d.addCallback(self._gotJob)
        d.addErrback(self._catchFailure)
        return d

    def getJobOutput(self):
        if self.remoteobj == None:
            return
        d = self.remoteobj.callRemote('getJobOutputByID', self._job.id)
        d.addCallback(self._gotJobOutput)
        d.addErrback(self._catchFailure)
        return d

    def _gotJobOutput(self, output):
        self.output = output
        self._job.output = output

    def getJobResult(self):
        if self.remoteobj == None:
            return
        d = self.remoteobj.callRemote('getJobResultByID', self._job.id)
        d.addCallback(self._gotJobResult)
        d.addErrback(self._catchFailure)
        return d

    def _gotJobResult(self, result):
        self.result = result
        self._job.result = result

    def syncJob(self):
        if self.remoteobj == None:
            if LOG_LEVEL > 2:
                print 'self.remoteobj is None'
            return
        if self.status == 'completed':
            if LOG_LEVEL > 2:
                print 'Stopping syncJob'
            if self.syncJob_task:
                if self.syncJob_task.running:
                    self.syncJob_task.stop()
            return

        try:
            d = self.remoteobj.callRemote('syncJob', self._job.id)
        except pb.DeadReferenceError:
            if self.syncJob_task:
                if self.syncJob_task.running:
                    self.syncJob_task.stop()
            return

        d.addCallback(self._gotJob)
        d.addErrback(self._catchFailure)

    def write_result(self, filename):
        result_file = open(filename, 'w')

        # skip the first element since that is not the actual result
        for line in self.result:
            line = str(line)
            result_file.write(line)
        result_file.close()

    def kill(self):
        r"""
        Kills the current job.

        """

        d = self.remoteobj.callRemote('killJob', self._job.id)
        d.addCallback(self._killedJob)
        d.addErrback(self._catchFailure)
        return d

    def _killedJob(self, jobID):
        if jobID:
            if LOG_LEVEL > 2:
                print str(jobID) + ' was successfully killed.\r'
