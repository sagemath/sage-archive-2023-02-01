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
import sqlite3

from twisted.python import log

import sage.dsage.database.sql_functions as sql_functions
from sage.dsage.misc.config import get_conf

class MonitorDatabase(object):
    """
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
     anonymous BOOL DEFAULT 0,
     last_connection timestamp
    )
    """

    def __init__(self, test=False):
        self.conf = get_conf(type='monitordb')
        self.tablename = 'monitors'
        if test:
            self.db_file = 'monitordb-test.db'
        else:
            self.db_file = self.conf['db_file']
            if not os.path.exists(self.db_file):
                dir, file = os.path.split(self.db_file)
                if not os.path.isdir(dir):
                    os.mkdir(dir)
        self.log_level = self.conf['log_level']
        self.log_file = self.conf['log_file']
        self.con = sqlite3.connect(self.db_file,
                    detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
        self.con.text_factory = sqlite3.OptimizedUnicode

        if sql_functions.table_exists(self.con, self.tablename) is None:
            sql_functions.create_table(self.con,
                                       self.tablename,
                                       self.CREATE_MONITOR_TABLE)
            self.con.commit()

    def _set_parameter(self, uuid, key, value):
        query = """UPDATE monitors
        SET %s=?
        WHERE uuid=?""" % (key)
        cur = self.con.cursor()
        cur.execute(query, (value, uuid))
        self.con.commit()

    def set_anonymous(self, uuid, anonymous=True):
        return self._set_parameter(uuid, 'anonymous', anonymous)

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
        query = """SELECT uuid, hostname, ip, anonymous, sage_version, os FROM monitors
        WHERE uuid=?"""
        cur = self.con.cursor()
        cur.execute(query, (uuid,))
        result = cur.fetchone()
        if result is None:
            return result
        columns = [desc[0] for desc in cur.description]
        monitor = dict(zip(columns, result))
        for k, v in monitor.iteritems():
            if k == 'anonymous':
                monitor[k] = bool(v)

        return monitor

    def get_monitor_list(self):
        """
        Returns a list of connected monitors.

        """

        query = """SELECT uuid, hostname, ip, anonymous, sage_version, os FROM monitors
        WHERE connected"""
        cur = self.con.cursor()
        cur.execute(query)
        result = cur.fetchall()
        columns = [desc[0] for desc in cur.description]
        monitors = [dict(zip(columns, monitor)) for monitor in result]
        for monitor in monitors:
            for k, v in monitor.iteritems():
                if k == 'anonymous':
                    monitor[k] = bool(v)

        return monitors

    def set_connected(self, uuid, connected=True):
        """
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
        """
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
        """
        Returns the number of workers.

        Parameters:
        connected -- bool
        busy -- bool

        """

        if connected and not busy:
            query = """SELECT workers FROM monitors WHERE connected AND NOT busy"""
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
        """
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
        """
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