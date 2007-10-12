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
import sqlite3

from twisted.python import log

from BTrees import OOBTree
import transaction

from sage.dsage.database.job import Job
import sage.dsage.database.sql_functions as sql_functions
from sage.dsage.misc.constants import DSAGE_DIR
from sage.dsage.misc.constants import DELIMITER

class JobDatabase(object):
    """
    Implementation of the job database.
    Common methods between the implementations should go here.

    """

    def __init__(self, db_file=os.path.join(DSAGE_DIR, 'db', 'dsage.db'),
                 job_failure_threshold=3,
                 log_file=os.path.join(DSAGE_DIR, 'server.log'),
                 log_level=0, test=False):
        self.test = test
        if self.test:
            self.db_file = 'test_jobdb.db'
            self.db_file = 'dsage_test.db'
            self.log_file = 'dsage_test.log'
            self.log_level = 5
            self.prune_in_days = 10
            self.job_failure_threshold = 3
        else:
            self.db_file = db_file
            if not os.path.exists(self.db_file):
                dir, file = os.path.split(self.db_file)
                if dir == '':
                  pass
                elif not os.path.isdir(dir):
                    os.mkdir(dir)
            self.job_failure_threshold = job_failure_threshold
            self.log_file = log_file
            self.log_level = log_level

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
     priority INTEGER DEFAULT 5,
     type TEXT,
     failures INTEGER DEFAULT 0,
     creation_time timestamp NOT NULL,
     update_time timestamp,
     start_time timestamp,
     finish_time timestamp,
     cpu_time REAL,
     wall_time REAL,
     verifiable BOOL,
     private BOOL DEFAULT 0,
     timeout INTEGER DEFAULT 600,
     killed BOOL DEFAULT 0
    );
    """

    def __init__(self, db_file=None, job_failure_threshold=3,
                 log_file=os.path.join(DSAGE_DIR, 'server.log'),
                 log_level=0, test=False):
        JobDatabase.__init__(self, db_file=db_file,
                             job_failure_threshold=job_failure_threshold,
                             log_file=log_file, log_level=log_level,
                             test=test)
        self.tablename = 'jobs'
        self.con = sqlite3.connect(
                  self.db_file,
                  isolation_level=None,
                  detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
        sql_functions.optimize_sqlite(self.con)
        # Don't use this, slow!
        # self.con.text_factory = sqlite3.OptimizedUnicode
        self.con.text_factory = str
        if not sql_functions.table_exists(self.con, self.tablename):
            sql_functions.create_table(self.con,
                                       self.tablename,
                                       self.CREATE_JOBS_TABLE)

    def _shutdown(self):
        self.con.commit()
        self.con.close()

    def get_job(self, anonymous=False):
        """
        Returns the first unprocessed job of the highest priority.

        """

        if anonymous:
            query = """SELECT * FROM jobs
                       WHERE status = 'new' AND killed = 0 ORDER BY priority
                       LIMIT 1
                    """
        else:
            query = """SELECT * FROM jobs
                       WHERE status = 'new' AND killed = 0 ORDER BY priority
                       LIMIT 1
                    """
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
        """
        Returns a jdict given a job id.

        """

        query = """SELECT
                job_id,
                name,
                status,
                output,
                code,
                result,
                killed,
                verifiable,
                monitor_id,
                cpu_time,
                wall_time,
                start_time,
                timeout,
                priority,
                failures
                FROM jobs WHERE job_id = ?"""

        cur = self.con.cursor()
        cur.execute(query, (job_id,))
        jtuple = cur.fetchone()
        jdict = self.create_jdict(jtuple, cur.description)
        del jtuple

        return jdict

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

        cur = self.con.cursor()
        query = """UPDATE jobs
                   SET %s=?
                   WHERE job_id=?
                """ % (key)
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
        jdict -- sage.dsage.database.Job.jdict

        """

        try:
            job_id = jdict['job_id']
        except KeyError, msg:
            job_id = None

        if job_id is None:
            job_id = self.random_string()
            if self.log_level > 3:
                log.msg('[JobDB] Creating a new job with id:', job_id)
            query = """INSERT INTO jobs
                    (job_id, status, creation_time) VALUES (?, ?, ?)"""
            cur = self.con.cursor()
            cur.execute(query, (job_id, 'new', datetime.datetime.now()))
            # self.con.commit()

        for k, v in jdict.iteritems():
            try:
                self._update_value(job_id, k, v)
            except (sqlite3.InterfaceError,
                    sqlite3.OperationalError,
                    sqlite3.IntegrityError), msg:
                if self.log_level > 3:
                    log.msg(DELIMITER)
                    log.msg('The following error is probably NOT BAD')
                    log.msg('key: %s, value: %s' % (k, v))
                    log.msg(msg)
                    log.msg(DELIMITER)
                continue

        return job_id

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
        try:
            jdict['data'] = str(jdict['data'])
        except Exception, msg:
            pass
        try:
            jdict['result'] = str(jdict['result'])
        except Exception, msg:
            pass

        return jdict

    def get_killed_jobs_list(self):
        """
        Returns a list of jobs which have been marked as killed.

        """

        query = """SELECT
                   job_id,
                   status,
                   killed,
                   verifiable,
                   monitor_id,
                   failures,
                   update_time
                   FROM jobs
                   WHERE killed = 1 AND status <> 'completed'"""

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
        """
        Sets the value of killed for a job.

        """

        return self._update_value(job_id, 'killed', killed)

    def get_active_jobs(self):
        """
        Returns a list of jdicts whose status is 'processing'

        """

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