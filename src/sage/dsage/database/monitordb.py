##############################################################################
#
#  DSAGE: Distributed SAGE
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
##############################################################################

import datetime
import os
import ConfigParser
import sqlite3 as sqlite

from twisted.python import log

import sage.dsage.database.sql_functions as sql_functions

class MonitorDatabase(object):
    r"""
    This table keeps track of workers.

    """

    CREATE_MONITOR_TABLE = """CREATE TABLE monitors
    (
     uuid text NOT NULL UNIQUE,
     hostname TEXT,
     ip TEXT,
     workers INTEGER,
     sage_version text,
     os text,
     kernel_version TEXT,
     cpus INTEGER,
     cpu_speed INTEGER,
     cpu_model TEXT,
     mem_total INTEGER,
     mem_free INTEGER,
     connected BOOL,
     busy BOOL,
     last_connection timestamp
    )
    """

    def __init__(self, test=False):
        self._getconf()
        self.tablename = 'monitors'
        if test:
            self.db_file = 'monitordb-test.db'
        else:
            self.db_file = self.DB_FILE
            if not os.path.exists(self.db_file):
                dir, file = os.path.split(self.db_file)
                if not os.path.isdir(dir):
                    os.mkdir(dir)

        self.con = sqlite.connect(self.db_file,
                    detect_types=sqlite.PARSE_DECLTYPES|sqlite.PARSE_COLNAMES)
        self.con.text_factory = str

        if sql_functions.table_exists(self.con, self.tablename) is None:
            sql_functions.create_table(self.con,
                                       self.tablename,
                                       self.CREATE_MONITOR_TABLE)
            self.con.commit()

    def _getconf(self):
        self.DSAGE_DIR = os.path.join(os.getenv('DOT_SAGE'), 'dsage')
        # Begin reading configuration
        try:
            conf_file = os.path.join(self.DSAGE_DIR, 'server.conf')
            config = ConfigParser.ConfigParser()
            config.read(conf_file)

            # TODO: This needs to be changed to use db_file
            self.DB_FILE = os.path.expanduser(config.get('db', 'db_file'))
            self.LOG_FILE = config.get('db_log', 'log_file')
            self.LOG_LEVEL = config.getint('db_log', 'log_level')
        except:
            print "Error reading '%s', run dsage.setup()" % conf_file
            raise
        # End reading configuration

    def add_monitor(self, host_info):
        query = """INSERT INTO monitors
        (uuid,
         hostname,
         ip,
         workers,
         sage_version,
         os,
         kernel_version,
         cpus,
         cpu_speed,
         cpu_model,
         mem_total,
         mem_free)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """

        uuid = host_info['uuid']
        hostname = host_info['hostname']
        ip = host_info['ip']
        workers = host_info['workers']
        sage_version = host_info['sage_version']
        os = host_info['os']
        kernel_version = host_info['kernel_version']
        cpus = host_info['cpus']
        cpu_speed = host_info['cpu_speed']
        cpu_model = host_info['cpu_model']
        mem_total = host_info['mem_total']
        mem_free = host_info['mem_free']

        cur = self.con.cursor()
        cur.execute(query, (uuid, hostname, ip, workers, sage_version, os,
                            kernel_version, cpus, cpu_speed, cpu_model,
                            mem_total, mem_free))
        self.con.commit()

    def get_monitor(self, uuid):
        query = """SELECT * FROM monitors WHERE uuid=?"""
        cur = self.con.cursor()
        cur.execute(query, (uuid,))

        return cur.fetchone()

    def get_monitor_list(self):
        r"""
        Returns a list of connected monitors.

        """
        query = """SELECT uuid FROM monitors WHERE connected"""

        cur = self.con.cursor()
        cur.execute(query)

        return cur.fetchall()

    def set_connected(self, uuid, connected=True):
        r"""
        Sets the connected status of a worker.

        Parameters:
        uuid -- string
        connected -- bool

        """

        cur = self.con.cursor()
        if connected:
            query = """UPDATE monitors SET connected=1, last_connection=?
            WHERE uuid=?"""
            cur.execute(query, (datetime.datetime.now(), uuid))
        else:
            query = """UPDATE monitors SET connected=0 WHERE uuid=?"""
            cur.execute(query, (uuid,))

        self.con.commit()

    def set_busy(self, uuid, busy):
        r"""
        Sets whether or not a worker is doing a job.

        """

        if busy:
            query = """UPDATE monitors SET busy=1 WHERE uuid=?"""
        else:
            query = """UPDATE monitors SET busy=0 WHERE uuid=?"""

        cur = self.con.cursor()
        cur.execute(query, (uuid,))
        self.con.commit()

    def get_worker_count(self, connected, busy):
        r"""
        Returns the number of workers.

        Parameters:
        connected -- bool
        busy -- bool
        """

        if connected and not busy:
            query = """SELECT workers FROM monitors
            WHERE connected AND NOT busy"""
        elif connected and busy:
            query = """SELECT workers FROM monitors WHERE connected AND busy"""
        elif connected:
            query = """SELECT workers FROM monitors WHERE connected"""
        else:
            query = "SELECT workers FROM monitors WHERE NOT connected"

        cur = self.con.cursor()
        cur.execute(query)
        result = cur.fetchall()

        return sum(w[0] for w in result)

    def get_cpu_speed(self, connected=True, busy=True):
        r"""
        Returns the aggregate cpu speed in Mhz.

        Parameters:
        connected -- bool

        """

        if connected:
            query = """SELECT cpu_speed, workers FROM monitors
            WHERE connected AND busy"""
        else:
            query = """SELECT cpu_speed, workers FROM monitors"""

        cur = self.con.cursor()
        cur.execute(query)

        result = cur.fetchall()

        cpu_speed = sum([s[0]*s[1] for s in result])

        return cpu_speed

    def get_cpu_count(self, connected=True):
        r"""
        Returns the number of cpus that are available.

        Parameters:
        connected -- bool

        """

        if connected:
            query = """SELECT workers, cpus FROM monitors WHERE connected"""
        else:
            query = """SELECT workers, cpus FROM monitors"""

        cur = self.con.cursor()
        cur.execute(query)

        result = cur.fetchall()

        cpu_count = sum(min(s[0:2]) for s in result)

        return cpu_count