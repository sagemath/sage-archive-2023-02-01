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
import datetime
import os
import ConfigParser
import random

from twisted.python import log

from ZODB import FileStorage, DB
from BTrees import OOBTree
import transaction

import sqlite3 as sqlite

from sage.dsage.database.job import Job

DSAGE_DIR = os.path.join(os.getenv('DOT_SAGE'), 'dsage')
# Begin reading configuration
try:
    conf_file = os.path.join(DSAGE_DIR, 'server.conf')
    config = ConfigParser.ConfigParser()
    config.read(conf_file)

    DB_FILE = os.path.expanduser(config.get('db', 'db_file'))
    PRUNE_IN_DAYS = config.getint('db', 'prune_in_days')
    STALE_IN_DAYS = config.getint('db', 'stale_in_days')
    FAILURE_THRESHHOLD = config.getint('db', 'failure_threshhold')
    LOG_FILE = config.get('db_log', 'log_file')
    LOG_LEVEL = config.getint('db_log', 'log_level')
except:
    print "Error reading '%s', please run dsage.setup() or fix manually"%conf_file
    sys.exit(-1)
# End reading configuration

class JobDatabaseZODB(object):
    r"""
    Implementation of DSage's database using ZODB.

    Parameters:
    db_file -- filename of the database to use
    test -- set to true for unittesting purposes

    """

    def __init__(self, test=False, read_only=False):
        if test:
            self.db_file = 'test_db.db'
        else:
            self.db_file = DB_FILE
            if not os.path.exists(self.db_file):
                dir, file = os.path.split(self.db_file)
                if not os.path.isdir(dir):
                    os.mkdir(dir)
        if read_only:
            self.storage = FileStorage.FileStorage(self.db_file, read_only=True)
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

        if self.getJobByID('GLOBALS') == None:
            self.__init_globals()

    def _shutdown(self):
        r"""
        Shuts down the database by closing all DB connections and
        storages.

        """
        transaction.commit()
        self.dbconn.close()
        self.db.close()

    def getJobID(self):
        r"""
        Increments jobID by one and returns the current jobID

        """

        gdict = self.__retrieveGlobals()
        jobNUM = gdict['next_job_num']
        gdict['next_job_num'] += 1
        jobID = self.random_string(10) + '_' + '%s' % jobNUM
        self.storeJob(gdict)
        if LOG_LEVEL > 1:
            log.msg('[DB] Incremented job num to: ', jobNUM)

        return jobID

    def random_string(self, length):
        r"""
        Returns a random string

        Parameters:
            length -- the length of the string

        """

        random.seed()
        letters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
        l = list(10 * letters)
        random.shuffle(l)
        s = ''.join(random.sample(l, length))

        return s

    def __retrieveGlobals(self):
        return self.getJobByID('GLOBALS')

    def __init_globals(self):
        gdict = dict(next_job_num = 1)
        self.storeJob(gdict)

    def getJob(self):
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

    def getJobByID(self, jobID):
        r"""
        Returns a job given a job id.

        Parameters:
        jobID -- the job id (int)

        """
        if not self.hasJob(jobID):
            return None
        else:
            return self.jobdb[jobID]

    def getJobsByAuthor(self, author, is_active, job_name=None):
        r"""
        Returns jobs created by author.

        Parameters:
        author -- the author name (str)
        is_active -- when set to true, only return active jobs (bool)
        job_name -- the job name to return (default=None)

        """

        if is_active:
            if job_name:
                jobs = [job for job in self.getActiveJobs()
                        if (job.author == author and job.name == job_name)]
            else:
                jobs = [job for job in self.getActiveJobs()
                        if job.author == author]

        else:
            if job_name:
                jobs = [job for job in self.getJobsList()
                        if (job.author == author and job.name == job_name)]
            else:
                jobs = [job for job in self.getJobsList()
                        if job.author == author]

        if LOG_LEVEL > 3:
            log.msg('[JobDatabaseZODB, getJobsByAuthor] ', jobs)
        return jobs

    def storeJob(self, job):
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
                # log.msg('[DB, storeJob] Setting jobID to GLOBALS')
                jobID = 'GLOBALS'
        elif isinstance(job, Job):
            jobID = job.id
            if not isinstance(jobID, str):
                jobID = self.getJobID()
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
        if LOG_LEVEL > 0:
            log.msg("[DB] Stored job %s" % (jobID))

        return jobID

    def removeJob(self, jobID):
        r"""
        Removes a job from the database.

        Parameters:
         jobID -- job id (int)

        """

        try:
            del self.jobdb[jobID]
            if LOG_LEVEL > 0:
                log.msg('[JobDB] Removed job %s from jobdb' % (jobID))
            transaction.commit()
            return jobID
        except:
            log.msg('[JobDB] Failure removing a job.')
            raise

    def getActiveJobs(self):
        r"""
        Returns a list containing active jobs, i.e.
        jobs that are currentl being processed.

        """

        return [job for job in self.getJobsList()
                if job.status == 'processing']

    def hasJob(self, jobID):
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

    def getJobsList(self):
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

    def getCompletedJobsList(self):
        r"""
        Returns an ordered list of all the completed
        jobs in the database.

        """

        return [job for job in self.getJobsList()
                if job.status == 'completed']

    def getIncompleteJobsList(self):
        r"""
        Returns an ordered list of jobs marked as incomplete.

        """

        return [job for job in self.getJobsList()
                if job.status == 'incomplete']

    def getKilledJobsList(self):
        r"""
        Returns a list of killed jobs.

        """
        return [job for job in self.getJobsList() if job.killed]

    def newJob(self, job):
        r"""
        Stores a new job in the database.

        Parameters:
            job -- a Job object

        """

        jobID = self.getJobID()
        job.id = jobID

        return self.storeJob(job)

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

    def cleanOldJobs(self):
        r"""
        Cleans out jobs that are older than PRUNE_IN_DAYS days.

        """

        log.msg('[DatabasePruner] Cleaning out old jobs...')
        jobs = self.jobdb.getJobsList()
        for job in jobs:
            delta =  datetime.datetime.now() - job.updated_time
            if delta > datetime.timedelta(PRUNE_IN_DAYS):
                self.jobdb.removeJob(job.id)
                log.msg('[DatabasePruner, cleanOldJobs] Deleted job ',
                        job.id)

    def cleanFailedJobs(self):
        r"""
        Cleans out jobs which are marked as having failed.

        Failure threshhold is set in FAILURE_THRESHHOLD, if a job
        exceeds this threshhold, we remove it from the job database.

        """

        jobs = self.jobdb.getJobsList()

        for job in jobs:
            if job.failures > FAILURE_THRESHHOLD:
                self.jobdb.removeJob(job.id)

    def prune(self):
        self.cleanOldJobs()
        self.cleanFailedJobs()
