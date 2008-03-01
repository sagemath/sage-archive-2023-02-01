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
import sqlite3
from sqlalchemy import *

from twisted.python import log
from zope.interface import Interface
from zope.interface import implements

import sage.dsage.database.sql_functions as sql_functions
from sage.dsage.database.job import Job, expand_job
from sage.dsage.misc.constants import DELIMITER, SERVER_LOG
from sage.dsage.misc.misc import random_string


class IJobDatabase(Interface):
    def get_job(self, authenticated=False):
        """
        Returns the first unprocessed job of the highest priority.

        """


    def set_job_uuid(self, job_id, uuid):
        """
        Sets the unique identifier for a job.

        """


    def get_job_uuid(self, job_id):
        """
        Returns the unique identifier for a job.

        """

    def get_all_jobs(self):
        """
        Returns all jobs in the database.

        """


    def get_job(self):
        """
        Returns a unprocessed job from the job queue.

        """


    def set_killed(self, job_id, killed=True):
        """
        Sets the killed status of a job.

        """


    def store_jdict(self, jdict):
        """
        Stores a job in the database.

        """


    def get_jobs_by_username(self, username, active=True):
        """
        Returns jobs belonging to a user.

        """


    def has_job(self, job_id):
        """
        Checks to see if the database has a particular job.

        """



class JobDatabaseSA(object):
    implements(IJobDatabase)

    def __init__(self, Session):
        self.sess = Session()
        self.failure_threshold = 5

    def _shutdown(self):
        self.sess.close()

    def get_job(self):
        q = self.sess.query(Job).order_by(Job.c.creation_time.asc()).limit(1)
        query = q.filter_by(status='new', killed=False)

        return query.first()

    def get_job_count(self):
        c = self.sess.query(Job).count()

        return c

    def get_job_by_id(self, job_id):
        return self.sess.query(Job).filter_by(job_id=job_id).first()

    def has_job(self, job_id):
        query = self.sess.query(Job).filter_by(job_id=job_id)

        return query.one() > 0

    def set_job_uuid(self, job_id, uuid):
        job = self.sess.query(Job).filter_by(job_id=job_id).first()
        job.uuid = uuid
        self.sess.save_or_update(job)
        self.sess.commit()

    def get_job_uuid(self, job_id):
        job = self.sess.query(Job).filter_by(job_id=job_id).first()

        return job.uuid

    def update_job(self, job):
        """
        Takes a job object and updates it in the database.

        """

        self.sess.save_or_update(job)
        self.sess.commit()

    def store_jdict(self, jdict):
        try:
            job_id = jdict['job_id']
            assert job_id != None
            job = self.sess.query(Job).filter_by(job_id=job_id).first()
            for k,v in jdict.iteritems():
                setattr(job, k, v)
        except:
            job_id = random_string(length=10)
            jdict['job_id'] = job_id
            job = expand_job(jdict)
        self.update_job(job)

        return job_id

    def get_job_range(self, start, end):
        q = self.sess.query(Job).order_by(Job.c.update_time.desc())
        j = q[start:end].all()

        return j

    def get_n_jobs(self, n):
        q = self.sess.query(Job)
        j = q.order_by(Job.c.creation_time.desc()).limit(n).all()

        return j

    def get_all_jobs(self):
        jobs = self.sess.query(Job).order_by(Job.c.creation_time.desc()).all()

        return jobs

    def get_active_jobs(self):
        jobs = self.sess.query(Job).filter_by(status='processing').all()

        return [job._reduce() for job in jobs]

    def get_jobs_by_username(self, username, status):
        q = self.sess.query(Job)
        q = q.filter_by(username=username).filter_by(status=status)

        return q.all()

    def get_killed_jobs_list(self):
        return self.sess.query(Job).filter_by(killed=True).all()

    def set_killed(self, job_id, killed=True):
        job = self.sess.query(Job).filter_by(job_id=job_id).first()
        job.killed = killed
        job.status = 'killed'
        self.sess.save_or_update(job)
        self.sess.commit()

        return job

class JobDatabaseSQLite(object):
    """
    Implementation of DSage's database using SQLite.

    Parameters:

    AUTHORS:
    Yi Qiang
    Alex Clemesha

    """

    implements(IJobDatabase)

    def __init__(self, db_conn, failure_threshold=3,
                 log_level=0, log_file=SERVER_LOG):
        self.con = db_conn
        self.failure_threshold = failure_threshold
        self.log_file = log_file
        self.log_level = log_level
        self.tablename = 'jobs'


    def _shutdown(self):
        self.con.commit()


    def get_job(self):
        """
        Returns the first unprocessed job of the highest priority.

        """

        query = """SELECT * FROM jobs
                   WHERE status = 'new' AND killed = 0
                   ORDER BY priority, creation_time
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
        query = """SELECT * from jobs
                   ORDER by update_time DESC
                """
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
                worker_info,
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


    def store_jdict(self, jdict):
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
            job_id = random_string(length=10)
            if self.log_level > 3:
                log.msg('[JobDB] Creating a new job with id:', job_id)
            query = """INSERT INTO jobs
                    (job_id, status, creation_time) VALUES (?, ?, ?)"""
            cur = self.con.cursor()
            cur.execute(query, (job_id, 'new', datetime.datetime.now()))
            # self.con.commit()

        for k, v in jdict.iteritems():
            if k == 'worker_info':
                v = str(v)
            try:
                sql_functions.update_value(self.con, 'jobs', 'job_id',
                                           job_id, k, v)
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

        # Convert worker_info from a str into a dict
        try:
            jdict['worker_info'] = eval(jdict['worker_info'])
        except TypeError:
            jdict['worker_info'] = {}

        # Convert buffer objects back to string
        try:
            jdict['data'] = str(jdict['data'])
        except KeyError, msg:
            pass
        try:
            jdict['result'] = str(jdict['result'])
        except KeyError, msg:
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
                   worker_info,
                   priority,
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

        sql_functions.update_value(self.con, 'jobs', 'job_id', job_id,
                                   'status', 'killed')
        return sql_functions.update_value(self.con, 'jobs', 'job_id',
                                          job_id, 'killed', killed)


    def get_active_jobs(self):
        """
        Returns a list of jdicts whose status is 'processing'

        """

        return self._get_jobs_by_parameter('status', 'processing')