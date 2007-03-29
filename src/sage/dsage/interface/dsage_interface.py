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
import random
import glob
import ConfigParser
import copy
import cPickle
import zlib
import threading
import time

from twisted.spread import pb
from twisted.internet import reactor, defer, error, task
from twisted.cred import credentials
from twisted.conch.ssh import keys

from sage.dsage.database.job import Job, expand_job
from sage.dsage.twisted.pb import PBClientFactory
from sage.dsage.twisted.misc import blocking_call_from_thread
from sage.dsage.errors.exceptions import NoJobException, NotConnectedException

class DSageThread(threading.Thread):
    def run(self):
        if not reactor.running:
            reactor.run(installSignalHandlers=0)

class DSage(object):
    r"""
    This object represents a connection to the distributed SAGE server.
    """

    def __init__(self, server=None, port=8081, username=None,
                 pubkey_file=None, privkey_file=None):

        # We will read the values in from the conf file first and let the
        # user override the values stored in the conf file by keyword
        # parameters

        self._getconf()

        if server is None:
            self.server = self.SERVER
        else:
            self.server = server
        if port is None:
            self.port = self.PORT
        else:
            self.port = port
        if username is None:
            self.username = self.USERNAME
        else:
            self.username = username
        if pubkey_file is None:
            self.pubkey_file = self.PUBKEY_FILE
        else:
            self.pubkey_file = pubkey_file
        if privkey_file is None:
            self.privkey_file = self.PRIVKEY_FILE
        else:
            self.privkey_file = privkey_file

        self.remoteobj = None
        self.result = None

        # public key authentication information
        self.pubkey_str = keys.getPublicKeyString(filename=self.pubkey_file)
        # try getting the private key object without a passphrase first
        try:
            self.priv_key = keys.getPrivateKeyObject(
                                filename=self.privkey_file)
        except keys.BadKeyError:
            passphrase = self._getpassphrase()
            self.priv_key = keys.getPrivateKeyObject(
                            filename=self.privkey_file,
                            passphrase=passphrase)
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

        self.jobs = []

        self.connect()

    def __str__(self):
        self.check_connected()
        self.info_str = 'Connected to DSAGE server at: ' \
                    + self.server + ':' + str(self.port)
        return self.info_str + '\r'

    def __call__(self, cmd, globals_=None, job_name=None):
        cmd = ['ans = %s\n' % (cmd),
               'print ans\n'
               "save(ans, 'ans')\n"
               "DSAGE_RESULT = 'ans.sobj'\n"]

        return self.eval(''.join(cmd), globals_=globals_, job_name=job_name)

    def __getstate__(self):
        d = copy.copy(self.__dict__)
        d['remoteobj'] = None
        return d

    def _getconf(self):
        # randomly generated string we will use to sign
        self.DATA =  ''.join([chr(i) for i in [random.randint(65, 123) for n in
                        range(500)]])
        self.DSAGE_DIR = os.path.join(os.getenv('DOT_SAGE'), 'dsage')
        # Begin reading configuration
        try:
            conf_file = os.path.join(self.DSAGE_DIR, 'client.conf')
            config = ConfigParser.ConfigParser()
            config.read(conf_file)

            self.LOG_FILE = config.get('log', 'log_file')
            self.LOG_LEVEL = config.getint('log', 'log_level')
            self.SSL = config.getint('ssl', 'ssl')
            self.USERNAME = config.get('auth', 'username')
            self.PRIVKEY_FILE = os.path.expanduser(config.get('auth',
                                                              'privkey_file'))
            self.PUBKEY_FILE = os.path.expanduser(config.get('auth',
                                                             'pubkey_file'))
            self.SERVER = config.get('general', 'server')
            self.PORT = config.getint('general', 'port')
        except Exception, msg:
            print msg
            raise

    def _getpassphrase(self):
        import getpass
        passphrase = getpass.getpass('Passphrase (Hit enter for None): ')

        return passphrase

    def _catch_failure(self, failure):
        if failure.check(error.ConnectionRefusedError):
            print 'Remote server refused the connection.'
            return
        print "Error: ", failure.getErrorMessage()
        print "Traceback: ", failure.printTraceback()

    def _connected(self, remoteobj):
        if self.LOG_LEVEL > 0:
            print 'Connected to remote server.\r'
        self.remoteobj = remoteobj
        self.remoteobj.notifyOnDisconnect(self._disconnected)

    def _disconnected(self, remoteobj):
        print 'Lost connection to the server.'
        self.info_str = 'Not connected.'

    def _got_my_jobs(self, jobs, job_name):
        if jobs == None:
            raise NoJobException
        if job_name:
            return [JobWrapper(self.remoteobj, job)
                    for job in jobs if job.name == job_name]

    def _killed_job(self, job_id):
        if job_id:
            if self.LOG_LEVEL > 2:
                print str(job_id) + ' was successfully killed.'

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
        factory = PBClientFactory()

        if self.SSL == 1:
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
                            self._catch_failure)

    def disconnect(self):
        print 'Disconnecting from server.'
        self.remoteobj = None

    def eval(self, cmd, globals_=None, job_name=None):
        r"""
        eval evaluates a command

        Parameters:
        cmd -- the sage command to be evaluated (str)
        globals -- a dict (see help for python's eval method)
        job_name -- an alphanumeric job name

        """

        self.check_connected()
        if not job_name or not isinstance(job_name, str):
            job_name = 'default job'

        type_ = 'sage'

        job = Job(id_=None, code=cmd, name=job_name,
                  user_id=self.username, type_=type_)

        wrapped_job = JobWrapper(self.remoteobj, job)
        if globals_ is not None:
            for k, v in globals_.iteritems():
                job.attach(k, v)

        return wrapped_job

    def eval_file(self, fname, job_name):
        r"""
        eval_file allows you to evaluate the contents of an entire file.

        Parameters:
            fname -- file name of the file you wish to evaluate

        """

        self.check_connected()

        type_ = 'file'
        cmd = open(fname).read()
        job = Job(id_=None, code=cmd, name=job_name,
                  user_id=self.username, type_=type_)

        wrapped_job = JobWrapper(self.remoteobj, job)

        return wrapped_job

    def send_job(self, job):
        r"""
        Sends a Job object to the server.

        """

        if not isinstance(job, Job):
            raise TypeError
        wrapped_job = JobWrapper(self.remoteobj, job)
        return wrapped_job

    def _got_job_id(self, id, job):
        job.id = id
        job.user_id = self.username

        self.jobs.append(job)

        pickled_job = job.pickle()
        d = self.remoteobj.callRemote('submit_job', pickled_job)
        d.addErrback(self._catch_failure)
        # d.addCallback(self._submitted, job)

        return JobWrapper(self.remoteobj, job)

    def eval_dir(self, dir, job_name):
        self.check_connected()
        os.chdir(dir)
        files = glob.glob('*.spyx')
        deferreds = []
        for file in files:
            sage_cmd = open(file).readlines()
            d = self.remoteobj.callRemote('get_next_job_id')
            d.addCallback(self._got_id, sage_cmd, job_name, file=True,
                          type_='spyx')
            d.addErrback(self._catch_failure)
            deferreds.append(d)
        d_list = defer.DeferredList(deferreds)
        return d_list

    def kill(self, job_id):
        r"""
        Kills a job given the job id.

        Parameters:
        job_id -- job id

        """

        d = self.remoteobj.callRemote('kill_job', job_id)
        d.addCallback(self._killed_job)
        d.addErrback(self._catch_failure)

    def get_my_jobs(self, is_active=False, job_name=None):
        r"""
        This method returns a list of jobs that belong to you.

        Parameters:
        is_active -- set to true to get only active jobs (bool)

        Use this method if you get disconnected from the server and wish to
        retrieve your old jobs back.

        """

        self.check_connected()

        d = self.remoteobj.callRemote('get_jobs_by_user_id',
                                      self.username,
                                      is_active,
                                      job_name)
        d.addCallback(self._got_my_jobs, job_name)
        d.addErrback(self._catch_failure)

        return d

    def cluster_speed(self):
        r"""
        Returns the speed of the cluster.

        """

        self.check_connected()

        return self.remoteobj.callRemote('get_cluster_speed')

    def check_connected(self):
        if self.remoteobj == None:
            raise NotConnectedException
        if self.remoteobj.broker.disconnected:
            raise NotConnectedException

class BlockingDSage(DSage):
    r"""This is the blocking version of DSage
    """

    def __init__(self, server=None, port=8081, username=None,
                 pubkey_file=None, privkey_file=None):
        # We will read the values in from the conf file first and let the
        # user override the values stored in the conf file by keyword
        # parameters

        self._getconf()

        if server is None:
            self.server = self.SERVER
        else:
            self.server = server
        if port is None:
            self.port = self.PORT
        else:
            self.port = port
        if username is None:
            self.username = self.USERNAME
        else:
            self.username = username
        if pubkey_file is None:
            self.pubkey_file = self.PUBKEY_FILE
        else:
            self.pubkey_file = pubkey_file
        if privkey_file is None:
            self.privkey_file = self.PRIVKEY_FILE
        else:
            self.privkey_file = privkey_file

        self.remoteobj = None
        self.result = None

        # public key authentication information
        self.pubkey_str = keys.getPublicKeyString(filename=self.pubkey_file)

        # try getting the private key object without a passphrase first
        try:
            self.priv_key = keys.getPrivateKeyObject(
                                filename=self.privkey_file)
        except keys.BadKeyError:
            passphrase = self._getpassphrase()
            self.priv_key = keys.getPrivateKeyObject(
                            filename=self.privkey_file,
                            passphrase=passphrase)

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

        self.jobs = []
        self.dsage_thread = DSageThread()
        self.dsage_thread.start()
        self.connect()

    def connect(self):
        r"""
        This methods establishes the conection to the remote server.

        """

        # TODO: Send a useful 'mind' object with the login request!
        # factory = pb.PBClientFactory()
        self.factory = PBClientFactory()

        if self.SSL == 1:
            from twisted.internet import ssl
            contextFactory = ssl.ClientContextFactory()
            blocking_call_from_thread(reactor.connectSSL,
                                   self.server, self.port,
                                   self.factory, contextFactory)
        else:
            blocking_call_from_thread(reactor.connectTCP,
                                   self.server, self.port,
                                   self.factory)

        d = self.factory.login(self.creds, None)
        d.addCallback(self._connected)
        d.addErrback(self._catch_failure)

        return d

    def eval(self, cmd, globals_=None, job_name=None, async=False):
        r"""
        eval evaluates a command

        Parameters:
        cmd -- the sage command to be evaluated (str)
        globals -- a dict (see help for python's eval method)
        job_name -- an alphanumeric job name
        async -- whether to use the async implementation of the method

        """

        self.check_connected()
        if not job_name or not isinstance(job_name, str):
            job_name = 'default_job'

        type_ = 'sage'

        job = Job(id_=None, code=cmd, name=job_name,
                  user_id=self.username, type_=type_)

        if globals_ is not None:
            for k, v in globals_.iteritems():
                job.attach(k, v)

        if async:
            wrapped_job = JobWrapper(self.remoteobj, job)
        else:
            wrapped_job = BlockingJobWrapper(self.remoteobj, job)

        return wrapped_job

    def send_job(self, job, async=False):
        r"""
        Sends a Job object to the server.

        Parameters:
        job -- a Job object to send to the remote server
        async -- if True, use async method of doing remote task

        """

        if not isinstance(job, Job):
            raise TypeError
        if async:
            wrapped_job = JobWrapper(self.remoteobj, job)
        else:
            wrapped_job = BlockingJobWrapper(self.remoteobj, job)

        return wrapped_job

    def get_my_jobs(self):
        r"""
        This method returns a list of jobs that belong to you.

        Parameters:
        is_active -- set to true to get only active jobs (bool)

        Use this method if you get disconnected from the server and wish to
        retrieve your old jobs back.

        """

        self.check_connected()

        jdicts = blocking_call_from_thread(self.remoteobj.callRemote,
                                           'get_jobs_by_user_id',
                                           self.username)

        return [expand_job(jdict) for jdict in jdicts]


    def cluster_speed(self):
        r"""
        Returns the speed of the cluster.

        """

        self.check_connected()

        return blocking_call_from_thread(self.remoteobj.callRemote,
                                      'get_cluster_speed')

    def list_monitors(self):
        r"""Returns a list of monitors connected to the server.

        """

        self.check_connected()

        return blocking_call_from_thread(self.remoteobj.callRemote,
                                         'get_monitor_list')

    def list_clients(self):
        r"""
        Returns a list of clients connected to the server.
        """

        self.check_connected()

        return blocking_call_from_thread(self.remoteobj.callRemote,
                                         'get_client_list')

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

        # d = self.remoteobj.callRemote('get_next_job_id')
        d = self.remoteobj.callRemote('submit_job', job.reduce())
        d.addCallback(self._got_jdict)
        d.addErrback(self._catch_failure)

    def __repr__(self):
        if self._job.status == 'completed' and not self._job.output:
            return 'No output.'
        elif not self._job.output:
            return 'No output yet.'

        return self._job.output

    def __getstate__(self):
        d = copy.copy(self.__dict__)
        d['remoteobj'] = None
        d['sync_job_task'] = None
        return d

    def _update_job(self, job):
        # This sets all the attributes of our JobWrapper object to match the
        # attributes of a Job object
        for k, v in Job.__dict__.iteritems():
            if isinstance(v, property):
                setattr(self, k, getattr(job, k))

        for k, v in job.__dict__.iteritems():
            setattr(self, k, getattr(job, k))

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

    def _catch_failure(self, failure):
        if failure.check(pb.DeadReferenceError, error.ConnectionLost):
            print 'Disconnected from server.'
        else:
            print "Error: ", failure.getErrorMessage()
            print "Traceback: ", failure.printTraceback()

    def _got_job(self, job):
        if job == None:
            return
        self._job = expand_job(job)
        self._update_job(self._job)

    def _got_jdict(self, jdict):
        self._job = expand_job(jdict)
        self.id = jdict['job_id']
        self._update_job(self._job)

    def get_job(self):
        if self.remoteobj is None:
            raise NotConnectedException
        if self.id is None:
            return
        d = self.remoteobj.callRemote('get_job_by_id', self.id)
        d.addCallback(self._got_job)
        d.addErrback(self._catch_failure)

        return d

    def get_job_output(self):
        if self.remoteobj == None:
            return
        d = self.remoteobj.callRemote('get_job_output_by_id', self._job.id)
        d.addCallback(self._got_job_output)
        d.addErrback(self._catch_failure)

        return d

    def _got_job_output(self, output):
        self.output = output
        self._job.output = output

    def get_job_result(self):
        if self.remoteobj == None:
            return
        d = self.remoteobj.callRemote('get_job_result_by_id', self._job.id)
        d.addCallback(self._got_job_result)
        d.addErrback(self._catch_failure)

        return d

    def _got_job_result(self, result):
        self.result = result
        self._job.result = result

    def sync_job(self):
        if self.remoteobj == None:
            if self.LOG_LEVEL > 2:
                print 'self.remoteobj is None'
            return
        if self.status == 'completed':
            if self.LOG_LEVEL > 2:
                print 'Stopping sync_job'
            if self.sync_job_task:
                if self.sync_job_task.running:
                    self.sync_job_task.stop()
            return

        try:
            d = self.remoteobj.callRemote('sync_job', self._job.id)
        except pb.DeadReferenceError:
            if self.sync_job_task:
                if self.sync_job_task.running:
                    self.sync_job_task.stop()
            return

        d.addCallback(self._got_job)
        d.addErrback(self._catch_failure)

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

        if self.id is not None:
            d = self.remoteobj.callRemote('kill_job', self.id)
            d.addCallback(self._killed_job)
            d.addErrback(self._catch_failure)
            return d
        else:
            return

    def _killed_job(self, job_id):
        return
        # if job_id:
        #     if self.LOG_LEVEL > 2:
        #         print str(job_id) + ' was successfully killed.\r'

class BlockingJobWrapper(JobWrapper):
    r"""
    Blocking version of the JobWrapper object.  This is to be used
    interactively.

    """

    def __init__(self, remoteobj, job):
        self.remoteobj = remoteobj
        self._job = job

        self._update_job(job)
        self.worker_info = self._job.worker_info

        jdict = blocking_call_from_thread(self.remoteobj.callRemote,
                                          'submit_job', job.reduce())
        self._job = expand_job(jdict)

    def __repr__(self):
        self.get_job()
        if self.status == 'completed' and not self.output:
            return 'No output.'
        elif not self.output:
            return 'No output yet.'

        return self.output

    def get_job(self):
        if self.remoteobj == None:
           raise NotConnectedException
        if self.status == 'completed':
            return

        job = blocking_call_from_thread(self.remoteobj.callRemote,
                                     'get_job_by_id', self._job.id)

        self._update_job(expand_job(job))

    def async_get_job(self):
        return JobWrapper.get_job(self)

    def kill(self):
        r"""
        Kills the current job.

        """

        job_id = blocking_call_from_thread(self.remoteobj.callRemote,
                                           'kill_job', self._job.id)
        return job_id

    def async_kill(self):
        r"""
        async version of kill

        """

        d = self.remoteobj.callRemote('kill_job', self.id)
        d.addCallback(self._killed_job)
        d.addErrback(self._catch_failure)
        return d

    def wait(self):
        timeout = 0.1
        while self.result is None:
            time.sleep(timeout)
            self.get_job()