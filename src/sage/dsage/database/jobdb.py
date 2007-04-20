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

import datetime
import os
import random
import string
import time
import sqlite3

from twisted.python import log

from ZODB import FileStorage, DB
from BTrees import OOBTree
import transaction

from sage.dsage.database.job import Job
import sage.dsage.database.sql_functions as sql_functions
from sage.dsage.misc.config import get_conf

class JobDatabase(object):
    """
    Implementation of the job database.
    Common methods between the implementations should go here.

    """

    def __init__(self):
        self.conf = get_conf(type='jobdb')
        self.db_file = self.conf['db_file']
        self.job_failure_threshold = int(self.conf['job_failure_threshold'])
        self.log_file = self.conf['log_file']
        self.log_level = int(self.conf['log_level'])
        self.prune_in_days = int(self.conf['prune_in_days'])

    def random_string(self, length=10):
        """
        Returns a random string

        Parameters:
            length -- the length of the string

        """

        random.seed()
        l = list((length*length) * (string.letters + string.digits))
        random.shuffle(l)
        s = ''.join(random.sample(l, length))

        return s

class JobDatabaseZODB(JobDatabase):
    """
    Implementation of DSage's database using ZODB.

    Parameters:
    db_file -- filename of the database to use
    test -- set to true for unittesting purposes

    """

    # The following effectively turns of the ZODB logger, which is OK for us.
    # Without this, one gets this annoying error message a lot:
    #       No handlers could be found for logger "ZODB.FileStorage"
    import logging
    logging.getLogger("ZODB.FileStorage").setLevel(10000000)
    logging.getLogger("ZODB.lock_file").setLevel(10000000)
    logging.getLogger("ZODB.Connection").setLevel(10000000)

    def __init__(self, test=False, read_only=False):
        JobDatabase.__init__(self)
        if test:
            self.db_file = 'test_db.db'
        else:
            self.db_file = self.DB_FILE
            if not os.path.exists(self.db_file):
                dir, file = os.path.split(self.db_file)
                if not os.path.isdir(dir):
                    os.mkdir(dir)
        if read_only:
            self.storage = FileStorage.FileStorage(self.db_file,
                                                   read_only=True)
        else:
            self.storage = FileStorage.FileStorage(self.db_file)
        self.db = DB(self.storage)
        self.dbconn = self.db.open()
        self.dbroot = self.dbconn.root()

        # main job database
        if not self.dbroot.has_key('jobdb'):
            self.dbroot['jobdb'] = OOBTree.OOBTree()
        # killed job database
        if not self.dbroot.has_key('killed_jobdb'):
            self.dbroot['killed_jobdb'] = OOBTree.OOBTree()
        # completed job database
        if not self.dbroot.has_key('completed_jobdb'):
            self.dbroot['completed_jobdb'] = OOBTree.OOBTree()
        self.jobdb = self.dbroot['jobdb']
        self.killed_jobdb = self.dbroot['killed_jobdb']
        self.completed_jobdb = self.dbroot['completed_jobdb']

        if self.get_job_by_id('GLOBALS') == None:
            self.__init_globals()

    def _shutdown(self):
        """
        Shuts down the database by closing all DB connections and
        storages.

        """
        transaction.commit()
        self.dbconn.close()
        self.db.close()

    def get_next_job_id(self):
        """
        Increments job_id by one and returns the current job_id

        """

        gdict = self.__retrieve_globals()
        jobNUM = gdict['next_job_num']
        gdict['next_job_num'] += 1
        job_id = self.random_string(10) + '_' + '%s' % jobNUM
        self.store_job(gdict)
        if self.log_level > 1:
            log.msg('[DB] Incremented job num to: ', jobNUM)

        return job_id

    def __retrieve_globals(self):
        return self.get_job_by_id('GLOBALS')

    def __init_globals(self):
        gdict = dict(next_job_num = 1)
        self.store_job(gdict)

    def get_job(self):
        """
        Returns the first non completed job in the job database.

        This method iterates through the database and returns the first job
        it finds which has its status set to None.  This does not guarantee
        that jobs are always processed in the right order.

        """
        for job_id, job in self.jobdb.iteritems():
            try:
                if isinstance(job, Job) and job.status == 'new':
                    return job
            except(KeyError):
                return None

    def get_job_by_id(self, job_id):
        """
        Returns a job given a job id.

        Parameters:
        job_id -- the job id (int)

        """
        if not self.has_job(job_id):
            return None
        else:
            return self.jobdb[job_id]

    def get_jobs_by_username(self, username, is_active, job_name=None):
        """
        Returns jobs created by author.

        Parameters:
        username -- the username name (str)
        is_active -- when set to true, only return active jobs (bool)
        job_name -- the job name to return (default=None)

        """

        if is_active:
            if job_name:
                jobs = [job for job in self.get_active_jobs()
                        if (job.username == username and job.name == job_name)]
            else:
                jobs = [job for job in self.get_active_jobs()
                        if job.username == username]

        else:
            if job_name:
                jobs = [job for job in self.get_jobs_list()
                        if (job.username == username and job.name == job_name)]
            else:
                jobs = [job for job in self.get_jobs_list()
                        if job.username == username]

        if self.log_level > 3:
            log.msg('[JobDatabaseZODB, get_jobs_by_username] ', jobs)
        return jobs

    def store_job(self, job):
        """
        Stores a job in the job database.

        Parameters:
        job -- the Job object to store (Job)

        Returns:
        the job id of the job stored (int)

        """

        # This is a hack to check to see if the job we are storing
        # is the GLOBAL dictionary which keeps track of our job ids
        if isinstance(job, dict):
            if job.has_key('next_job_num'):
                # log.msg('[DB, store_job] Setting job_id to GLOBALS')
                job_id = 'GLOBALS'
        elif isinstance(job, Job):
            job_id = job.job_id
            if not isinstance(job_id, str):
                job_id = self.get_next_job_id()
            job.update_time = datetime.datetime.now()
        else:
            raise TypeError

        self.jobdb[job_id] = job
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        # Need this hack to notify ZODB that the dict was modified
        # WITHOUT IT CHANGES WILL NOT BE WRITTEN BACK TO DISK
        if job_id != 'GLOBALS':
            job._p_changed = True
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        transaction.commit()
        if self.log_level > 0:
            log.msg("[DB] Stored job %s" % (job_id))

        return job_id

    def remove_job(self, job_id):
        """
        Removes a job from the database.

        Parameters:
         job_id -- job id (int)

        """

        try:
            del self.jobdb[job_id]
            if self.log_level > 0:
                log.msg('[JobDB] Removed job %s from jobdb' % (job_id))
            transaction.commit()
            return job_id
        except:
            log.msg('[JobDB] Failure removing a job.')
            raise

    def get_active_jobs(self):
        """
        Returns a list containing active jobs, i.e.
        jobs that are currentl being processed.

        """

        return [job for job in self.get_jobs_list()
                if job.status == 'processing']

    def has_job(self, job_id):
        """
        Checks if the database contains a job corresponding to job id.

        Parameters:
        job_id -- the job id (int)

        """

        # For some reason has_key returns 0/1 versus True/False
        b = self.jobdb.has_key(job_id)
        if b > 0:
            return True
        else:
            return False

    def get_jobs_list(self):
        """
        Returns an ordered list of all the jobs in
        the database.

        """

        sorted_list = []

        for job_id, job in self.jobdb.iteritems():
            if job_id == 'GLOBALS': # Skip the GLOBALS job
                continue
            sorted_list.append((job.job_id, job))
        sorted_list.sort()

        return [job[1] for job in sorted_list]

    def get_completed_jobs_list(self):
        """
        Returns an ordered list of all the completed
        jobs in the database.

        """

        return [job for job in self.get_jobs_list()
                if job.status == 'completed']

    def get_in_complete_jobs_list(self):
        """
        Returns an ordered list of jobs marked as incomplete.

        """

        return [job for job in self.get_jobs_list()
                if job.status == 'incomplete']

    def get_killed_jobs_list(self):
        """
        Returns a list of killed jobs.

        """
        return [job for job in self.get_jobs_list() if job.killed]

    def new_job(self, job):
        """
        Stores a new job in the database.

        Parameters:
            job -- a Job object

        """

        job_id = self.get_next_job_id()
        job.job_id = job_id

        return self.store_job(job)

class JobDatabaseSQLite(JobDatabase):
    """
    Implementation of DSage's database using SQLite.

    Parameters:
    test -- set to true for unittesting purposes

    AUTHORS:
    Yi Qiang
    Alex Clemesha

    """

    # TODO: SQLite does *NOT* enforce foreign key constraints
    # Must do manual checking.

    CREATE_JOBS_TABLE = """CREATE TABLE jobs
    (job_id TEXT NOT NULL UNIQUE,
     name TEXT,
     username TEXT REFERENCES clients(username),
     monitor_id TEXT REFERENCES monitors(uuid),
     worker_info TEXT,
     code TEXT,
     data BLOB,
     output TEXT,
     result BLOB,
     status TEXT NOT NULL,
     priority INTEGER DEFAULT 10,
     type TEXT,
     failures INTEGER DEFAULT 0,
     creation_time timestamp NOT NULL,
     update_time timestamp,
     finish_time timestamp,
     verifiable BOOL,
     timeout INTEGER DEFAULT 600,
     killed BOOL DEFAULT 0
    );
    """

    def __init__(self, test=False):
        JobDatabase.__init__(self)
        if test:
            self.db_file = 'test_jobdb.db'
        else:
            if not os.path.exists(self.db_file):
                dir, file = os.path.split(self.db_file)
                if not os.path.isdir(dir):
                    os.mkdir(dir)

        self.tablename = 'jobs'
        self.con = sqlite3.connect(
                   self.db_file,
                   detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
        self.con.text_factory = sqlite3.OptimizedUnicode # Don't want unicode objects
        if not sql_functions.table_exists(self.con, self.tablename):
            sql_functions.create_table(self.con,
                                       self.tablename,
                                       self.CREATE_JOBS_TABLE)

    def __add_test_data(self):
        INSERT_JOB = """INSERT INTO jobs
        (uid,
         username,
         code,
         data,
         output,
         worker_id,
         status,
         priority,
         creation_time,
         update_time,
         finish_time
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
        """

        dn = lambda: datetime.datetime.now() # short hand
        D = [(self.random_string(), 1, None, None, None, 1,
             'new', 1, dn(), dn(), None),
             (self.random_string(), 1, None, None, None, 1,
             0, 1, dn(), dn(), None),
             (self.random_string(), 1, None, None, None, 1,
             0, 1, dn(), dn(), None)]

        log.msg('Sleeping for 0.5 second to test timestamps...')
        time.sleep(0.5)

        D.extend([(self.random_string(), 1, None, None, None, 1,
                   0, 1, dn(), dn(), None),
                  (self.random_string(), 1, None, None, None, 1,
                   0, 1, dn(), dn(), None),
                  (self.random_string(), 1, None, None, None, 1,
                   0, 1, dn(), dn(), None)])

        self.con.executemany(INSERT_JOB, D)
        self.con.commit()

    def _shutdown(self):
        self.con.commit()
        self.con.close()

    def get_job(self, anonymous=False):
        """
        Returns the first unprocessed job of the highest priority.

        """

        if anonymous:
            query = "SELECT * FROM jobs WHERE status = 'new' AND killed = 0"
        else:
            query = "SELECT * FROM jobs WHERE status = 'new' AND killed = 0"
        cur = self.con.cursor()
        cur.execute(query)
        jtuple = cur.fetchone()

        return self.create_jdict(jtuple, cur.description)

    def set_job_uuid(self, job_id, uuid):
        query = "UPDATE jobs SET monitor_id=? WHERE job_id=?"
        cur = self.con.cursor()
        cur.execute(query, (uuid, job_id))
        self.con.commit()

    def get_all_jobs(self):
        query = "SELECT * from jobs"
        cur = self.con.cursor()
        cur.execute(query)
        result = cur.fetchall()

        return [self.create_jdict(jtuple, cur.description)
                for jtuple in result]

    def get_job_by_id(self, job_id):
        query = "SELECT * FROM jobs WHERE job_id = ?"
        cur = self.con.cursor()
        cur.execute(query, (job_id,)) # Need to cast it to int for SAGE
        jtuple = cur.fetchone()
        return self.create_jdict(jtuple, cur.description)

    def _get_jobs_by_parameter(self, key, value):
        """
        Returns a particular result given a key and value.

        """

        query = """SELECT * FROM jobs where %s = ?""" % (key)
        cur = self.con.cursor()
        cur.execute(query, (value,))
        result = cur.fetchall()

        return [self.create_jdict(jtuple, cur.description)
                for jtuple in result]

    def _set_parameter(self, job_id, key, value):
        query = """UPDATE jobs
        SET %s=?
        WHERE job_id=?""" % (key)
        cur = self.con.cursor()
        cur.execute(query, (value, job_id))
        self.con.commit()

    def _update_value(self, job_id, key, value):
        """
        Sets the appropriate value for a job in the database.

        """

        query = """UPDATE jobs
        SET %s=?
        WHERE job_id=?
        """ % (key)
        cur = self.con.cursor()
        if key == 'data' or key == 'result': # Binary objects
            if value != None:
                cur.execute(query, (sqlite3.Binary(value), job_id))
        else:
            cur.execute(query, (value, job_id))
        self.con.commit()

    def store_job(self, jdict):
        """
        Stores a job based on information from Job.jdict.
        The keys of the dictionary should correspond to the columns in the
        'jobs' table.

        Parameters:
        job -- a sage.dsage.database.Job object

        """

        job_id = jdict['job_id']

        if job_id is None:
            job_id = self.random_string()
            if self.log_level > 3:
                log.msg('[JobDB] Creating a new job with id:', job_id)
            query = """INSERT INTO jobs
                    (job_id, status, creation_time) VALUES (?, ?, ?)"""
            cur = self.con.cursor()
            cur.execute(query, (job_id, 'new', datetime.datetime.now()))
            self.con.commit()

        for k, v in jdict.iteritems():
            try:
                self._update_value(job_id, k, v)
            except (sqlite3.InterfaceError,
                    sqlite3.OperationalError,
                    sqlite3.IntegrityError), msg:
                if self.log_level > 3:
                    log.msg('key: %s, value: %s' % (k, v))
                    log.msg(msg)
                continue

        return self.get_job_by_id(job_id)

    def create_jdict(self, jtuple, row_description):
        """
        Creates a jdict out of a job_tuple.

        """

        if jtuple is None:
            return None

        columns = [desc[0] for desc in row_description]
        jdict = dict(zip(columns, jtuple))

        # Convert 0/1 into python booleans
        jdict['killed'] = bool(jdict['killed'])
        jdict['verifiable'] = bool(jdict['verifiable'])

        # Convert buffer objects back to string
        jdict['data'] = str(jdict['data'])
        jdict['result'] = str(jdict['result'])

        return jdict

    def get_killed_jobs_list(self):
        """
        Returns a list of jobs which have been marked as killed.

        """
        query = "SELECT * from jobs where killed = 1 AND status <> 'completed'"
        cur = self.con.cursor()
        cur.execute(query)
        killed_jobs = cur.fetchall()
        return [self.create_jdict(jdict, cur.description)
                for jdict in killed_jobs]

    def get_jobs_by_username(self, username, active=True):
        """
        Returns a list of jobs belonging to 'username'

        """

        if active:
            query = """SELECT * from jobs WHERE username = ?
            AND status IN ('new', 'processing')"""
        else:
            query = """SELECT * from jobs WHERE username = ? """
        cur = self.con.cursor()
        cur.execute(query, (username,))

        return [self.create_jdict(jtuple, cur.description) for jtuple in cur]

    def has_job(self, job_id):
        """
        Checks if the database contains a job with the given uid.

        """

        job = self.get_job_by_id(job_id)
        if job is None:
            return False
        else:
            return True

    def set_killed(self, job_id, killed=True):
        self._update_value(job_id, 'killed', killed)

    def get_active_jobs(self):
        return self._get_jobs_by_parameter('status', 'processing')

class DatabasePruner(object):
    """
    DatabasePruner is responsible for cleaning out the database.

    """
    def __init__(self, jobdb):
        """
        Parameters:
            jobdb -- a JobDatabase object

        """

        self.jobdb = jobdb

    def clean_old_jobs(self):
        """
        Cleans out jobs that are older than PRUNE_IN_DAYS days.

        """

        log.msg('[DatabasePruner] Cleaning out old jobs...')
        jobs = self.jobdb.get_jobs_list()
        for job in jobs:
            delta =  datetime.datetime.now() - job.update_time
            if delta > datetime.timedelta(self.jobdb.prune_in_days):
                self.jobdb.remove_job(job.job_id)
                log.msg('[DatabasePruner, clean_old_jobs] Deleted job ',
                        job.job_id)

    def clean_failed_jobs(self):
        """
        Cleans out jobs which are marked as having failed.

        """

        jobs = self.jobdb.get_jobs_list()

        for job in jobs:
            if job.failures > self.jobdb.job_failure_threshold:
                self.jobdb.remove_job(job.job_id)

    def prune(self):
        self.clean_old_jobs()
        self.clean_failed_jobs()
