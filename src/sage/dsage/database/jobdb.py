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
import ConfigParser
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

class JobDatabase(object):
    r"""
    Implementation of the job database.
    Common methods between the implementations should go here.

    """

    def __init__(self):
        self._getconf()

    def _getconf(self):
        self.DSAGE_DIR = os.path.join(os.getenv('DOT_SAGE'), 'dsage')
        # Begin reading configuration
        try:
            conf_file = os.path.join(self.DSAGE_DIR, 'server.conf')
            config = ConfigParser.ConfigParser()
            config.read(conf_file)

            self.DB_FILE = os.path.expanduser(config.get('db', 'db_file'))
            self.PRUNE_IN_DAYS = config.getint('db', 'prune_in_days')
            self.STALE_IN_DAYS = config.getint('db', 'stale_in_days')
            self.JOB_FAILURE_THRESHOLD = config.getint('db',
                                                       'job_failure_threshold')
            self.LOG_FILE = config.get('db_log', 'log_file')
            self.LOG_LEVEL = config.getint('db_log', 'log_level')
        except:
            print "Error reading '%s', run dsage.setup()" % conf_file
            raise
        # End reading configuration

    def random_string(self, length=10):
        r"""
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
    r"""
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
        r"""
        Shuts down the database by closing all DB connections and
        storages.

        """
        transaction.commit()
        self.dbconn.close()
        self.db.close()

    def get_next_job_id(self):
        r"""
        Increments jobID by one and returns the current jobID

        """

        gdict = self.__retrieve_globals()
        jobNUM = gdict['next_job_num']
        gdict['next_job_num'] += 1
        jobID = self.random_string(10) + '_' + '%s' % jobNUM
        self.store_job(gdict)
        if self.LOG_LEVEL > 1:
            log.msg('[DB] Incremented job num to: ', jobNUM)

        return jobID

    def __retrieve_globals(self):
        return self.get_job_by_id('GLOBALS')

    def __init_globals(self):
        gdict = dict(next_job_num = 1)
        self.store_job(gdict)

    def get_job(self):
        r"""
        Returns the first non completed job in the job database.

        This method iterates through the database and returns the first job
        it finds which has its status set to None.  This does not guarantee
        that jobs are always processed in the right order.

        """
        for jobID, job in self.jobdb.iteritems():
            try:
                if isinstance(job, Job) and job.status == 'new':
                    return job
            except(KeyError):
                return None

    def get_job_by_id(self, jobID):
        r"""
        Returns a job given a job id.

        Parameters:
        jobID -- the job id (int)

        """
        if not self.has_job(jobID):
            return None
        else:
            return self.jobdb[jobID]

    def get_jobs_by_user_id(self, user_id, is_active, job_name=None):
        r"""
        Returns jobs created by author.

        Parameters:
        user_id -- the user_id name (str)
        is_active -- when set to true, only return active jobs (bool)
        job_name -- the job name to return (default=None)

        """

        if is_active:
            if job_name:
                jobs = [job for job in self.get_active_jobs()
                        if (job.user_id == user_id and job.name == job_name)]
            else:
                jobs = [job for job in self.get_active_jobs()
                        if job.user_id == user_id]

        else:
            if job_name:
                jobs = [job for job in self.get_jobs_list()
                        if (job.user_id == user_id and job.name == job_name)]
            else:
                jobs = [job for job in self.get_jobs_list()
                        if job.user_id == user_id]

        if self.LOG_LEVEL > 3:
            log.msg('[JobDatabaseZODB, get_jobs_by_user_id] ', jobs)
        return jobs

    def store_job(self, job):
        r"""
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
                # log.msg('[DB, store_job] Setting jobID to GLOBALS')
                jobID = 'GLOBALS'
        elif isinstance(job, Job):
            jobID = job.id
            if not isinstance(jobID, str):
                jobID = self.get_next_job_id()
            job.update_time = datetime.datetime.now()
        else:
            raise TypeError

        self.jobdb[jobID] = job
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        # Need this hack to notify ZODB that the dict was modified
        # WITHOUT IT CHANGES WILL NOT BE WRITTEN BACK TO DISK
        if jobID != 'GLOBALS':
            job._p_changed = True
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        transaction.commit()
        if self.LOG_LEVEL > 0:
            log.msg("[DB] Stored job %s" % (jobID))

        return jobID

    def remove_job(self, jobID):
        r"""
        Removes a job from the database.

        Parameters:
         jobID -- job id (int)

        """

        try:
            del self.jobdb[jobID]
            if self.LOG_LEVEL > 0:
                log.msg('[JobDB] Removed job %s from jobdb' % (jobID))
            transaction.commit()
            return jobID
        except:
            log.msg('[JobDB] Failure removing a job.')
            raise

    def get_active_jobs(self):
        r"""
        Returns a list containing active jobs, i.e.
        jobs that are currentl being processed.

        """

        return [job for job in self.get_jobs_list()
                if job.status == 'processing']

    def has_job(self, jobID):
        r"""
        Checks if the database contains a job corresponding to job id.

        Parameters:
        jobID -- the job id (int)

        """

        # For some reason has_key returns 0/1 versus True/False
        b = self.jobdb.has_key(jobID)
        if b > 0:
            return True
        else:
            return False

    def get_jobs_list(self):
        r"""
        Returns an ordered list of all the jobs in
        the database.

        """

        sorted_list = []

        for jobID, job in self.jobdb.iteritems():
            if jobID == 'GLOBALS': # Skip the GLOBALS job
                continue
            sorted_list.append((job.num, job))
        sorted_list.sort()

        return [job[1] for job in sorted_list]

    def get_completed_jobs_list(self):
        r"""
        Returns an ordered list of all the completed
        jobs in the database.

        """

        return [job for job in self.get_jobs_list()
                if job.status == 'completed']

    def get_in_complete_jobs_list(self):
        r"""
        Returns an ordered list of jobs marked as incomplete.

        """

        return [job for job in self.get_jobs_list()
                if job.status == 'incomplete']

    def get_killed_jobs_list(self):
        r"""
        Returns a list of killed jobs.

        """
        return [job for job in self.get_jobs_list() if job.killed]

    def new_job(self, job):
        r"""
        Stores a new job in the database.

        Parameters:
            job -- a Job object

        """

        jobID = self.get_next_job_id()
        job.id = jobID

        return self.store_job(job)

class JobDatabaseSQLite(JobDatabase):
    r"""
    Implementation of DSage's database using SQLite.

    Parameters:
    test -- set to true for unittesting purposes

    """

    # TODO: SQLite does *NOT* enforce foreign key constraints
    # Must do manual checking.

    CREATE_JOBS_TABLE = """CREATE TABLE jobs
    (id TEXT NOT NULL UNIQUE,
     name TEXT,
     user_id INTEGER REFERENCES users(id),
     worker_id INTEGER REFERENCES workers(id),
     worker_info TEXT,
     code TEXT,
     data BLOB,
     output TEXT,
     result BLOB,
     status TEXT NOT NULL,
     priority TEXT,
     type TEXT,
     failures INTEGER DEFAULT 0,
     creation_time timestamp NOT NULL,
     update_time timestamp,
     finish_time timestamp,
     verifiable BOOL,
     killed BOOL DEFAULT 0
    );
    """

    def __init__(self, test=False):
        JobDatabase.__init__(self)
        if test:
            self.db_file = 'test_jobdb.db'
        else:
            self.db_file = self.DB_FILE
            if not os.path.exists(self.db_file):
                dir, file = os.path.split(self.db_file)
                if not os.path.isdir(dir):
                    os.mkdir(dir)

        self.tablename = 'jobs'
        self.con = sqlite3.connect(
                   self.db_file,
                   detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
        self.con.text_factory = str # Don't want unicode objects
        if not sql_functions.table_exists(self.con, self.tablename):
            sql_functions.create_table(self.con,
                                       self.tablename,
                                       self.CREATE_JOBS_TABLE)

    def __add_test_data(self):
        INSERT_JOB = """INSERT INTO jobs
        (uid,
         user_id,
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

        print 'Sleeping for 0.5 second to test timestamps...'
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
        r"""
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

    def get_all_jobs(self):
        query = "SELECT * from jobs"
        cur = self.con.cursor()
        cur.execute(query)
        return cur.fetchall()

    def get_job_by_id(self, job_id):
        query = "SELECT * FROM jobs WHERE id = ?"
        cur = self.con.cursor()
        cur.execute(query, (job_id,)) # Need to cast it to int for SAGE
        jtuple = cur.fetchone()
        return self.create_jdict(jtuple, cur.description)

    def get_job_by_keywords(self, **kwargs):
        raise NotImplementedError

    def update_value(self, id, key, value):
        r"""
        Sets the appropriate value for a job in the database.

        """

        query = """UPDATE jobs
        SET %s=?
        WHERE id=?
        """ % (key)
        cur = self.con.cursor()
        if key == 'data' or key == 'result': # Binary objects
            if value is not None:
                cur.execute(query, (sqlite3.Binary(value), id))
        else:
            cur.execute(query, (value, id))
        self.con.commit()

    def store_job(self, jdict):
        r"""
        Stores a job based on information from Job.jdict.
        The keys of the dictionary should correspond to the columns in the
        'jobs' table.

        Parameters:
        job -- a sage.dsage.database.Job object

        """

        id = jdict['id']

        if id is None:
            id = self.random_string()
            if self.LOG_LEVEL > 3:
                log.msg('[JobDB] Creating a new job with id:', id)
            query = """INSERT INTO jobs
                    (id, status, creation_time) VALUES (?, ?, ?)"""
            cur = self.con.cursor()
            cur.execute(query, (id, 'new', datetime.datetime.now()))
            self.con.commit()

        for k, v in jdict.iteritems():
            try:
                self.update_value(id, k, v)
            except (sqlite3.InterfaceError,
                    sqlite3.OperationalError,
                    sqlite3.IntegrityError), msg:
                if self.LOG_LEVEL > 3:
                    print k, v
                    print msg
                continue

        return self.get_job_by_id(id)

    def create_jdict(self, jtuple, row_description):
        r"""
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
        r"""
        Returns a list of jobs which have been marked as killed.

        """
        query = "SELECT * from jobs where killed = 1 AND status <> 'completed'"
        cur = self.con.cursor()
        cur.execute(query)
        killed_jobs = cur.fetchall()
        return [self.create_jdict(jdict, cur.description)
                for jdict in killed_jobs]

    def get_jobs_by_user_id(self, user_id):
        r"""
        Returns a list of jobs belonging to 'user_id'

        """

        query = """SELECT * from jobs where user_id = ? """
        cur = self.con.cursor()
        cur.execute(query, (user_id,))

        return [self.create_jdict(jtuple, cur.description) for jtuple in cur]

    def has_job(self, id):
        r"""
        Checks if the database contains a job with the given uid.

        """

        job = self.get_job_by_id(id)
        if job is None:
            return False
        else:
            return True

class DatabasePruner(object):
    r"""
    DatabasePruner is responsible for cleaning out the database.

    """
    def __init__(self, jobdb):
        r"""
        Parameters:
            jobdb -- a JobDatabase object

        """

        self.jobdb = jobdb

    def clean_old_jobs(self):
        r"""
        Cleans out jobs that are older than PRUNE_IN_DAYS days.

        """

        log.msg('[DatabasePruner] Cleaning out old jobs...')
        jobs = self.jobdb.get_jobs_list()
        for job in jobs:
            delta =  datetime.datetime.now() - job.update_time
            if delta > datetime.timedelta(self.jobdb.PRUNE_IN_DAYS):
                self.jobdb.remove_job(job.id)
                log.msg('[DatabasePruner, clean_old_jobs] Deleted job ',
                        job.id)

    def clean_failed_jobs(self):
        r"""
        Cleans out jobs which are marked as having failed.

        """

        jobs = self.jobdb.get_jobs_list()

        for job in jobs:
            if job.failures > self.jobdb.JOB_FAILURE_THRESHOLD:
                self.jobdb.remove_job(job.id)

    def prune(self):
        self.clean_old_jobs()
        self.clean_failed_jobs()
