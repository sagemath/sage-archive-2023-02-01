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

import os
import glob
import copy
import cPickle
import zlib
import threading
import time
from getpass import getuser

from twisted.cred.credentials import Anonymous
from twisted.internet.threads import blockingCallFromThread
from twisted.internet import reactor

from sage.dsage.database.job import Job, expand_job
from sage.dsage.misc.misc import random_str
from sage.dsage.misc.constants import DSAGE_DIR

class DSageThread(threading.Thread):
    """
    DSage thread

    """

    def run(self):
        if not reactor.running:
            try:
                reactor.run(installSignalHandlers=0)
            except AttributeError, msg:
                pass
                # This is a temporary workaround for a weird bug in reactor
                # during shutdown that one sees doing doctests (on some
                # systems?).


class DSage(object):
    """
    This object represents a connection to the distributed SAGE server.
    Parameters:
    server -- str
    port -- int
    username -- str
    pubkey_file -- str (Default: ~/.sage/dsage/dsage_key.pub)
    privkey_file -- str (Default: ~/.sage/dsage/dsage_key)
    log_level -- int (Default: 0)
    ssl -- int (Default: 1)
    """

    def __init__(self, server='localhost', port=8081,
                 username=getuser(),
                 pubkey_file=os.path.join(DSAGE_DIR, 'dsage_key.pub'),
                 privkey_file=os.path.join(DSAGE_DIR, 'dsage_key'),
                 log_level=0, ssl=True, testing=False):

        from twisted.cred import credentials
        from twisted.conch.ssh import keys
        from twisted.spread import banana
        banana.SIZE_LIMIT = 100*1024*1024 # 100 MegaBytes

        self.server = server
        self.port = port
        self.username = username
        self._data = random_str(500)
        self.ssl = ssl
        self._log_level = log_level
        self._pubkey_file = pubkey_file
        self._privkey_file = privkey_file
        self._remoteobj = None
        self.result = None
        self.info_str = 'Connected to: %s:%s'
        self._testing = testing

        if not self._testing:
            self._pubkey = keys.Key.fromFile(self._pubkey_file)
            try:
                self._privkey = keys.Key.fromFile(self._privkey_file)
            except keys.BadKeyError, msg:
                passphrase = self._getpassphrase()
                self._privkey = keys.Key.fromFile(self._privkey_file,
                                                  passphrase=passphrase)
            self._algorithm = 'rsa'
            self._blob = self._pubkey.blob()
            self._signature = self._privkey.sign(self._data)
            self._creds = credentials.SSHPrivateKey(self.username,
                                                   self._algorithm,
                                                   self._blob,
                                                   self._data,
                                                   self._signature)
        else:
            self.username = 'tester'
        self.connect()


    def __repr__(self):
        return self.__str__()

    def __str__(self):
        if self.is_connected():
            return self.info_str % (self.server, self.port)
        else:
            return 'Not connected.'

    def __call__(self, cmd, user_vars=None, load_files=[], job_name=None):
        cmd = ['ans = %s\n' % (cmd),
               'print ans\n',
               "DSAGE_RESULT = ans\n"]

        return self.eval(''.join(cmd), user_vars=user_vars,
                                       load_files=load_files,
                                       job_name=job_name)

    def __getstate__(self):
        d = copy.copy(self.__dict__)
        d['remoteobj'] = None

        return d

    def _getpassphrase(self):
        import getpass
        passphrase = getpass.getpass('Passphrase (Hit enter for None): ')

        return passphrase

    def _catch_failure(self, failure):
        print "Error connecting: %s" % failure.getErrorMessage()

    def _connected(self, remoteobj):
        if self._log_level > 0:
            print 'Connected to remote server.\r'
        self._remoteobj = remoteobj
        self._remoteobj.notifyOnDisconnect(self._disconnected)

    def _disconnected(self, remoteobj):
        print '[DSage] Closed connection to %s' % (self.server)
        self.info_str = 'Not connected.'

    def _got_my_jobs(self, jobs, job_name):
        from sage.dsage.errors.exceptions import NoJobException
        if jobs == None:
            raise NoJobException
        if job_name:
            return [JobWrapper(self._remoteobj, job)
                    for job in jobs if job.name == job_name]

    def _killed_job(self, job_id):
        pass

    def restore(self, remoteobj):
        """
        This method restores a connection to the server.

        """

        self._remoteobj = remoteobj

    def connect(self):
        """
        This methods establishes the conection to the remote server.

        """

        from twisted.internet import reactor
        from sage.dsage.twisted.pb import ClientFactory
        factory = ClientFactory(self._login, (), {})
        factory.continueTrying = False # Do not attempt to reconnect

        if self.ssl == 1:
            # Old, uses OpenSSL, SAGE uses GNUTLS now
            # from twisted.internet import ssl
            # contextFactory = ssl.ClientContextFactory()
            # reactor.connectSSL(self.server,
            #                    self.port,
            #                    factory,
            #                    contextFactory)
            from gnutls.interfaces.twisted import X509Credentials
            cred = X509Credentials()
            reactor.connectTLS(self.server, self.port, factory, cred)
        else:
            reactor.connectTCP(self.server, self.port, factory)

    def _login(self, *args, **kwargs):
        if self._testing:
            d = self.factory.login(Anonymous(), None)
        else:
            d = self.factory.login(self._creds, None)
        d.addCallback(self._connected)
        d.addErrback(self._catch_failure)

        return d

    def disconnect(self):
        print 'Disconnecting from server.'
        self._remoteobj = None

    def eval(self, cmd, timeout=0, user_vars=None, job_name=None):
        """
        eval evaluates a command

        Parameters:
        cmd -- the sage command to be evaluated (str)
        globals -- a dict (see help for python's eval method)
        job_name -- an alphanumeric job name

        """

        self.is_connected()
        if not job_name or not isinstance(job_name, str):
            job_name = 'default job'

        kind = 'sage'

        # We have to convert timeout to a python int so it will not cause
        # security exceptions with twisted.

        job = Job(job_id=None, code=cmd, name=job_name,
                  username=self.username, timeout=timeout, kind=kind)

        wrapped_job = JobWrapper(self._remoteobj, job)
        if user_vars is not None:
            for k, v in user_vars.iteritems():
                job.attach(k, v)

        return wrapped_job

    def eval_file(self, fname, job_name, async=False):
        """
        eval_file allows you to evaluate the contents of an entire file.

        Parameters:
            fname -- file name of the file you wish to evaluate

        """

        self.is_connected()

        kind = 'file'
        cmd = open(fname).read()
        job = Job(job_id=None, code=cmd, name=job_name,
                  username=self.username, kind=kind)

        if async:
            wrapped_job = JobWrapper(self._remoteobj, job)
        else:
            wrapped_job = BlockingJobWrapper(self._remoteobj, job)

        return wrapped_job

    def send_job(self, job):
        """
        Sends a Job object to the server.

        """

        if not isinstance(job, Job):
            raise TypeError
        wrapped_job = JobWrapper(self._remoteobj, job)
        return wrapped_job

    def _got_job_id(self, id_, job):
        job.job_id = id_
        job.username = self.username
        pickled_job = job.pickle()
        d = self._remoteobj.callRemote('submit_job', pickled_job)
        d.addErrback(self._catch_failure)
        # d.addCallback(self._submitted, job)

        return JobWrapper(self._remoteobj, job)

    def eval_dir(self, dir_, job_name):
        from twisted.internet import defer
        self.is_connected()
        os.chdir(dir_)
        files = glob.glob('*.spyx')
        deferreds = []
        for file_ in files:
            sage_cmd = open(file_).readlines()
            d = self._remoteobj.callRemote('get_next_job_id')
            d.addCallback(self._got_job_id, sage_cmd, job_name, file=True,
                          kind='spyx')
            d.addErrback(self._catch_failure)
            deferreds.append(d)
        d_list = defer.DeferredList(deferreds)
        return d_list

    def kill(self, job_id, async=False):
        """
        Kills a job given the job id.

        Parameters:
        job_id -- job id

        """

        if async:
            d = self._remoteobj.callRemote('kill_job', job_id)
            d.addCallback(self._killed_job)
            d.addErrback(self._catch_failure)
        else:
            job_id = blockingCallFromThread(self._remoteobj.callRemote,
                                               'kill_job',
                                               job_id)


    def get_my_jobs(self, is_active=False, job_name=None):
        """
        This method returns a list of jobs that belong to you.

        Parameters:
        is_active -- set to true to get only active jobs (bool)

        Use this method if you get disconnected from the server and wish to
        retrieve your old jobs back.

        """

        self.is_connected()

        d = self._remoteobj.callRemote('get_jobs_by_username',
                                      self.username,
                                      is_active,
                                      job_name)
        d.addCallback(self._got_my_jobs, job_name)
        d.addErrback(self._catch_failure)

        return d

    def cluster_speed(self):
        """
        Returns the speed of the cluster.

        """

        self.is_connected()

        return self._remoteobj.callRemote('get_cluster_speed')

    def is_connected(self):
        if self._remoteobj == None:
            return False
        if self._remoteobj.broker.disconnected:
            raise False
        return True

class BlockingDSage(DSage):
    """
    This is the blocking version of the DSage interface.

    """
    def __init__(self, server='localhost', port=8081,
                 username=getuser(),
                 pubkey_file=os.path.join(DSAGE_DIR, 'dsage_key.pub'),
                 privkey_file=os.path.join(DSAGE_DIR, 'dsage_key'),
                 log_level=0, ssl=True, testing=False):
        self._dsage_thread = DSageThread()
        self._dsage_thread.setDaemon(False)
        self._dsage_thread.start()
        DSage.__init__(self, server=server, port=port, username=username,
                       pubkey_file=pubkey_file, privkey_file=privkey_file,
                       log_level=log_level, ssl=ssl, testing=testing)


    def connect(self):
        """
        This methods establishes the conection to the remote server.

        """

        from twisted.internet import reactor
        from sage.dsage.twisted.pb import ClientFactory

        self.factory = ClientFactory(self._login, (), {})
        self.factory.continueTrying = False

        if self.ssl:
            from gnutls.interfaces.twisted import X509Credentials
            cred = X509Credentials()
            blockingCallFromThread(reactor, reactor.connectTLS, self.server,
                                   self.port, self.factory, cred)
        else:
            blockingCallFromThread(reactor, reactor.connectTCP, self.server, self.port,
                                   self.factory)

    def _login(self, *args, **kwargs):
        if self._testing:
            d = self.factory.login(Anonymous(), None)
        else:
            d = self.factory.login(self._creds, None)
        d.addCallback(self._connected)
        d.addErrback(self._catch_failure)

        return d

    def eval(self, cmd, user_vars=None, job_name=None, timeout=600,
             load_files=[], priority=5, async=False):
        """
        eval evaluates a command

        Parameters:
        cmd -- the sage command to be evaluated (str)
        user_vars -- a dict of predefined variables you want to use.
        job_name -- an alphanumeric job name
        timeout -- an upper limit on how long the job runs before the worker
                   restarts itself
        load_files -- list of files to load before executing the job
        priority -- priority of the job created (0-5)
        async -- whether to use the async implementation of the method

        """

        self.is_connected()
        kind = 'sage'

        job = Job(job_id=None, code=cmd, name=job_name,
                  username=self.username, timeout=timeout, priority=priority,
                  kind=kind)

        for fname in load_files:
            if os.path.exists(fname):
                job.attach_file(fname)

        if user_vars is not None:
            for k, v in user_vars.iteritems():
                job.attach(k, v)

        if async:
            wrapped_job = JobWrapper(self._remoteobj, job)
        else:
            wrapped_job = BlockingJobWrapper(self._remoteobj, job)

        return wrapped_job

    def send_job(self, job, async=False):
        """
        Sends a Job object to the server.

        Parameters:
        job -- a Job object to send to the remote server
        async -- if True, use async method of doing remote task

        """

        if not isinstance(job, Job):
            raise TypeError
        if async:
            wrapped_job = JobWrapper(self._remoteobj, job)
        else:
            wrapped_job = BlockingJobWrapper(self._remoteobj, job)

        return wrapped_job

    def get_my_jobs(self, status='new'):
        """
        This method returns a list of jobs that belong to you.

        Parameters:
        active -- set to true to get only active jobs (bool)

        Use this method if you get disconnected from the server and wish to
        retrieve your old jobs back.

        """

        self.is_connected()
        jdicts = blockingCallFromThread(reactor, self._remoteobj.callRemote,
                                        'get_jobs_by_username',
                                        self.username, status)

        return [expand_job(jdict) for jdict in jdicts]


    def kill_all(self):
        """
        Kills all of your active jobs.

        """

        active_jobs = self.get_my_jobs(active=True)

        for job in active_jobs:
            self.kill(job.job_id)

    def cluster_speed(self):
        """
        Returns the speed of the cluster.

        """

        self.is_connected()

        return blockingCallFromThread(reactor, self._remoteobj.callRemote,
                                         'get_cluster_speed')

    def get_workers_list(self):
        """Returns a list of monitors connected to the server.

        """

        self.is_connected()

        return blockingCallFromThread(reactor, self._remoteobj.callRemote,
                                         'get_worker_list')

    def get_client_list(self):
        """
        Returns a list of clients connected to the server.
        """

        self.is_connected()

        return blockingCallFromThread(reactor, self._remoteobj.callRemote,
                                         'get_client_list')

    def get_worker_count(self):
        """
        Returns the number of busy and free workers.

        """

        self.is_connected()

        return blockingCallFromThread(reactor, self._remoteobj.callRemote,
                                         'get_worker_count')

    def web_server_url(self):
        """
        Returns the web server url.
        """

        self.is_connected()

        return blockingCallFromThread(reactor, self._remoteobj.callRemote,
                                      'web_server_url')

    def web_view(self):
        """
        Opens the dsage server's web interface in a browser.

        """

        from sage.server.misc import open_page
        url = self.web_server_url()
        address = url.split(':')[1].strip('/')
        port = int(url.split(':')[2].strip('/'))
        open_page(address, port, False)

    def server_log(self, n=50):
        return blockingCallFromThread(reactor, self._remoteobj.callRemote,
                                      'read_log', n, 'server')

    def worker_log(self, n=50):
        return blockingCallFromThread(reactor, self._remoteobj.callRemote,
                                      'read_log', n, 'worker')


class JobWrapper(object):
    """
    Represents a remote job.

    Parameters:
        remoteobj -- the PB server's remoteobj
        job -- a Job object (job)

    """

    def __init__(self, remoteobj, job):
        self._remoteobj = remoteobj
        self._update_job(job._reduce())
        # d = self._remoteobj.callRemote('get_next_job_id')
        try:
            d = self._remoteobj.callRemote('submit_job', job._reduce())
        except Exception, msg:
            print msg
        d.addCallback(self._got_job_id)
        d.addCallback(self._got_jdict)
        d.addErrback(self._catch_failure)

    def __repr__(self):
        return self.job_id

    def __str__(self):
        if self.status == 'completed' and not self.output:
            return 'No output. (Done)'
        elif not self.output:
            return 'No output yet. (Not done)'

        return self.output

    def __getstate__(self):
        d = copy.copy(self.__dict__)
        d['remoteobj'] = None
        d['sync_job_task'] = None

        return d

    def _update_job(self, jdict):
        self._jdict = jdict
        job = expand_job(jdict)
        for k, v in job.__dict__.iteritems():
            setattr(self, k, v)

    def wait(self):
        from twisted.internet import reactor
        timeout = 0.5
        while self._job.result is None:
            reactor.iterate(timeout)

    def save(self, filename=None):
        if filename is None:
            filename = str(self._job.name)
        filename += '.sobj'
        f = open(filename, 'w')
        cPickle.dump(self, f, 2)

        return filename

    def restore(self, dsage):
        self._remoteobj = dsage.remoteobj

    def _catch_failure(self, failure):
        from twisted.internet import error
        from twisted.spread import pb
        if failure.check(pb.DeadReferenceError, error.ConnectionLost):
            print 'Disconnected from server.'
        else:
            pass
            # print "Error: ", failure.getErrorMessage()
            # print "Traceback: ", failure.printTraceback()

    def _got_job_id(self, job_id):
        self.job_id = job_id
        try:
            d = self._remoteobj.callRemote('get_job_by_id', job_id)
        except Exception, msg:
            raise

        return d

    def _got_jdict(self, jdict):
        self.job_id = jdict['job_id']
        self._update_job(jdict)

    def get_job(self):
        from sage.dsage.errors.exceptions import NotConnectedException

        if self._remoteobj is None:
            raise NotConnectedException
        if self.job_id is None:
            return
        try:
            d = self._remoteobj.callRemote('get_job_by_id', self.job_id)
        except Exception, msg:
            raise

        d.addCallback(self._got_jdict)
        d.addErrback(self._catch_failure)

        return d

    def get_job_output(self):
        if self._remoteobj == None:
            return
        try:
            d = self._remoteobj.callRemote('get_job_output_by_id',
                                           self.job_id)
        except Exception, msg:
            raise

        d.addCallback(self._got_job_output)
        d.addErrback(self._catch_failure)

        return d

    def _got_job_output(self, output):
        self.output = output

    def get_job_result(self):
        if self._remoteobj == None:
            return
        try:
            d = self._remoteobj.callRemote('get_job_result_by_id',
                                           self.job_id)
        except Exception, msg:
            raise

        d.addCallback(self._got_job_result)
        d.addErrback(self._catch_failure)

        return d

    def _got_job_result(self, result):
        self.result = result

    def sync_job(self):
        from twisted.spread import pb
        if self._remoteobj == None:
            # if self._log_level > 2:
            #     print 'self._remoteobj is None'
            return
        if self.status == 'completed':
            # if self._log_level > 2:
            #     print 'Stopping sync_job'
            if self.sync_job_task:
                if self.sync_job_task.running:
                    self.sync_job_task.stop()
            return

        try:
            d = self._remoteobj.callRemote('sync_job', self.job_id)
        except pb.DeadReferenceError:
            if self.sync_job_task:
                if self.sync_job_task.running:
                    self.sync_job_task.stop()
            return

        d.addCallback(self._got_jdict)
        d.addErrback(self._catch_failure)

    def write_result(self, filename):
        result_file = open(filename, 'w')

        # skip the first element since that is not the actual result
        for line in self.result:
            line = str(line)
            result_file.write(line)
        result_file.close()

    def kill(self):
        """
        Kills the current job.

        """

        if self.job_id is not None:
            try:
                d = self._remoteobj.callRemote('kill_job', self.job_id)
            except Exception, msg:
                print 'Unable to kill %s because %s'  % (self.job_id, msg)
                return
            d.addCallback(self._killed_job)
            d.addErrback(self._catch_failure)
            return d
        else:
            return

    def _killed_job(self, job_id):
        return

class BlockingJobWrapper(JobWrapper):
    """
    Blocking version of the JobWrapper object.  This is to be used
    interactively.

    """

    def __init__(self, remoteobj, job):
        self._update_job(job._reduce())
        self._remoteobj = remoteobj
        self.job_id = blockingCallFromThread(reactor, self._remoteobj.callRemote,
                                           'submit_job', job._reduce())

    def __repr__(self):
        if self.killed:
            return 'Job %s was killed' % (self.job_id)
        if self.status != 'completed':
            self.get_job()
        if self.status == 'completed' and not self.output:
            return 'No output.'
        if not self.output:
            return 'No output yet.'
        else:
            return self.output

    def get_job(self):
        from sage.dsage.errors.exceptions import NotConnectedException

        if self._remoteobj == None:
            raise NotConnectedException
        if self.status == 'completed':
            return

        jdict = blockingCallFromThread(reactor, self._remoteobj.callRemote,
                                        'get_job_by_id', self.job_id)

        self._update_job(jdict)

    def async_get_job(self):
        return JobWrapper.get_job(self)

    def rerun(self):
        """
        Resubmits the current job.

        """
        self.job_id = blockingCallFromThread(reactor, self._remoteobj.callRemote,
                                             'submit_job', self._jdict)
    def kill(self):
        """
        Kills the current job.

        """

        job_id = blockingCallFromThread(reactor, self._remoteobj.callRemote,
                                           'kill_job', self.job_id)
        self.job_id = job_id
        self.killed = True

        return job_id


    def async_kill(self):
        """
        async version of kill

        """

        d = self._remoteobj.callRemote('kill_job', self.job_id)
        d.addCallback(self._killed_job)
        d.addErrback(self._catch_failure)

        return d


    def wait(self, timeout=None):
        """
        Waits on a job until it is completed.

        Parameters:
        timeout -- number of seconds to wait, if it has not completed by then
                   it will raise RunTimeError if it is set to None,
                   it will wait indefinitely until the job is completed
        """

        import signal

        if timeout is None:
            while self.status != 'completed':
                # print 'Wating...'
                time.sleep(1.0)
                self.get_job()
        else:
            def handler(signum, frame):
                raise RuntimeError('Maximum wait time exceeded.')
            signal.signal(signal.SIGALRM, handler)
            signal.alarm(timeout)
            while self.status != 'completed':
                time.sleep(1.0)
                self.get_job()
            signal.alarm(0)
