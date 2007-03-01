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
import sqlite3 as sqlite

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

    def get_jobs_by_author(self, author, is_active, job_name=None):
        r"""
        Returns jobs created by author.

        Parameters:
        author -- the author name (str)
        is_active -- when set to true, only return active jobs (bool)
        job_name -- the job name to return (default=None)

        """

        if is_active:
            if job_name:
                jobs = [job for job in self.get_active_jobs()
                        if (job.author == author and job.name == job_name)]
            else:
                jobs = [job for job in self.get_active_jobs()
                        if job.author == author]

        else:
            if job_name:
                jobs = [job for job in self.get_jobs_list()
                        if (job.author == author and job.name == job_name)]
            else:
                jobs = [job for job in self.get_jobs_list()
                        if job.author == author]

        if self.LOG_LEVEL > 3:
            log.msg('[JobDatabaseZODB, get_jobs_by_author] ', jobs)
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
            job.updated_time = datetime.datetime.now()
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

    Author: Alex Clemesha (clemesha@gmail.com)

    Parameters:
    db_file -- filename of the database to use
    test -- set to true for unittesting purposes

    """

    # Queries and other often reused vars go here
    JOB_TABLE_COLUMNS = 12 #change this if more columns are added
    JOB_ROW_STRING = "(" + "?,"*JOB_TABLE_COLUMNS + ")"

    COLUMNS = ("num",
               "id",
               "file",
               "author",
               "data",
               "output",
               "worker_info",
               "status",
               "creation_time",
               "update_time",
               "finish_time")

    CREATE_NEW_TABLE = """CREATE TABLE jobs
    (num INTEGER PRIMARY KEY,
     id text NOT NULL UNIQUE,
     file text,
     author text NOT NULL,
     data text,
     output text,
     worker_info text,
     status text NOT NULL,
     creation_time timestamp NOT NULL,
     update_time timestamp NOT NULL,
     finish_time timestamp NOT NULL
    );
    """

    INSERT_JOB = """INSERT INTO jobs
    (num,
     id,
     file,
     author,
     data,
     output,
     worker_info,
     status,
     creation_time,
     update_time,
     finish_time
    )
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
    """

    def __init__(self, db_file=None, test=False):
        JobDatabase.__init__(self)
        if test:
            self.db_file = 'test_db.db'
        elif db_file:
            self.db_file = db_file
        else:
            self.db_file = 'jobdb.db'

        self.tablename = 'jobs'
        self.con = sqlite.connect(
                   self.db_file,
                   detect_types=sqlite.PARSE_DECLTYPES|sqlite.PARSE_COLNAMES)
        self.con.text_factory = str
        if not sql_functions.table_exists(self.con, self.tablename):
            sql_functions.create_table(self. con,
                                       self.tablename,
                                       self.CREATE_NEW_TABLE)

    def _shutdown(self):
        self.con.commit()
        self.con.close()

    def get_job_by_id(self, jobID, tablename=None):
        if tablename is None:
            tablename = self.tablename
        query = "SELECT * FROM %s WHERE id = '%s'" % (tablename, jobID)
        self.con.execute(query)
        jobtuple = self.con.fetchall()[0]
        return self.create_job(jobtuple)

    def get_job_by_num(self, num):
        query = "SELECT * FROM ? WHERE num = ?"
        self.con.execute(query, (self.tablename, num))
        job_tuple = self.con.fetchall()[0]
        # return self.create_job(jobtuple)
        return job_tuple

    # def create_job(self, jobtuple):
    #     job = Job()
    #     print jobtuple
    #     for k,v in zip(self.COLUMNS, jobtuple):
    #         print k,v
    #         print type(v)
    #         setattr(job, k, v)
    #     return job

    def store_job(self, job_dict):
        r"""
        Stores a job based on information from a dictionary.
        The keys of the dictionary should correspond to the columns in the
        'jobs' table.

        Parameters:
        job_dict -- dictionary

        """

        # TODO: Build a query from the job dictionary and store it in the
        # database

    def make_job_tuple(self):
        dn = lambda: datetime.datetime.now()
        D = (None, 'a', 'aaa', 'alex', dn(), dn(), dn(), dn())
        return D

    def add_test_data(self):
        dn = lambda: datetime.datetime.now()
        D = [(None, self.random_string(), 'aaa', 'alex',
              None, None, None, 'new', dn(), dn(), dn()),
             (None, self.random_string(), 'aaa', 'alex',
              None, None, None, 'new', dn(), dn(), dn()),
             (None, self.random_string(), 'aaa', 'alex',
              None, None, None, 'new', dn(), dn(), dn())]
        time.sleep(1)
        D.extend([(None, self.random_string(), 'aaa', 'alex',
                   None, None, None, 'new', dn(), dn(), dn()),
                  (None, self.random_string(), 'aaa', 'alex',
                   None, None, None, 'new', dn(), dn(), dn()),
                  (None, self.random_string(), 'aaa', 'alex',
                   None, None, None, 'new', dn(), dn(), dn())])

        query = self.INSERT_JOB % (self.tablename)
        self.con.executemany(query, D)
        self.con.commit()

    def insert_new_row_from_dict(self, tbl, dvalues):
        """
        """

        query = "INSERT INTO jobs VALUES (?, ?, ?, ?, ?, ?, ?, ?)"
        self.con.execute(query, dvalues)
        self.con.commit()

    def _general_select(self, table, condition, rowlist):
        """Contsruct a general SQL query.

        'SELECT *rows FROM table WHERE condition'

        Is this function a good idea?? Maybe more error checking??
        """
        if not isinstance(rowlist, (list, tuple)):
            return TypeError, "rowlist must be a list of a tuple"
        rl = [s.replace(' ', '') for s in rowlist]
        rl = ', '.join(rl)
        query =  "SELECT %s FROM %s WHERE %s"%(rl, str(table), str(condition))
        self.con.execute(query)
        return self.con.fetchall()

    def table_exists(self, table):
        """Check if a given table exists.

        If the below query is not None, then the table exists
        """
        query = """SELECT name FROM sqlite_master WHERE type = 'table' AND name = '%s' """ % (table)
        self.con.execute(query)
        cur = self.con.cursor()
        result = cur.fetchone()
        print "\n======= The status of '%s' is '%s' ========\n" % (table, result)
        return result

    def select_user_jobs(self, user):
        query = "SELECT num, creation_time FROM jobs WHERE user = '%s'"%str(user)
        self.con.execute(query)
        print "\nQUERY: '%s' "%query
        print "QUERY RESULTS: "
        for j, t in self.con.fetchall():
            print "job string: %s, timestamp: %s" % (str(j), str(t))

    def select_job_number(self, num):
        query = "SELECT num, creation_time FROM jobs WHERE %s > 2 " % str(num)
        self.con.execute(query)
        print "\nQUERY: '%s' "%query
        print "QUERY RESULTS: "
        for j, t in self.con.fetchall():
            print "job string: %s, timestamp: %s" % (str(j), str(t))

    def _newer_time_diff_seconds(self, older, newer, delta):
        """time difference
        find all time that are no older than 'delta' seconds.

        """
        d = datetime.timedelta(seconds=delta)
        diff = newer - d
        return older < diff

    def _get_oldest_job_datetime(self):
        query = "SELECT creation_time FROM jobs where num = 1"
        self.con.execute(query)
        return self.con.fetchone()[0]

    def select_job_date(self):
        OLDEST = self._get_oldest_job_datetime()
        query = "SELECT num, user, creation_time as 'creation_time [timestamp]' FROM jobs"
        self.con.execute(query)
        print "\nQUERY: '%s' "%query
        print "QUERY RESULTS: "
        for n, u, t in self.con.fetchall():
            if self._newer_time_diff_seconds(OLDEST, t, 2):
                print " %s, %s, %s"%(str(n), str(u), str(t))

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
            delta =  datetime.datetime.now() - job.updated_time
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
