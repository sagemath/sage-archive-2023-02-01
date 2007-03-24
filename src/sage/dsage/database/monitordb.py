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
import cStringIO
import xml.dom.minidom

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
            pass
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

    def get_cpu_speed(self, connected=True):
        r"""
        Returns the aggregate cpu speed in Mhz.

        Parameters:
        connected -- bool

        """

        if connected:
            query = """SELECT cpu_speed, workers FROM monitors WHERE
            connected"""
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
            query = """SELECT cpus FROM monitors WHERE connected"""
        else:
            query = """SELECT cpus FROM monitors"""

        cur = self.con.cursor()
        cur.execute(query)

        result = cur.fetchall()

        cpu_count = sum(s[0] for s in result)

        return cpu_count

    def generate_widget_xml(self):
        r"""
        This method returns a an XML document to be consumed by the Dashboard
        widget

        """

        def create_gauge(doc):
            gauge = doc.createElement('gauge')
            doc.appendChild(gauge)

            return doc, gauge

        def add_totalAgentCount(doc, gauge):
            totalAgentCount = doc.createElement('totalAgentCount')
            gauge.appendChild(totalAgentCount)
            working_workers = self.get_worker_count(connected=True, busy=True)
            free_workers = self.get_worker_count(connected=True, busy=False)
            disconnected_workers = self.get_worker_count(connected=False,
                                                         busy=False)
            total_workers = (working_workers +
                             free_workers +
                             disconnected_workers)
            count = doc.createTextNode(str(total_workers))
            totalAgentCount.appendChild(count)

            return doc, totalAgentCount

        def add_onlineAgentCount(doc, gauge):
            onlineAgentCount = doc.createElement('onlineAgentCount')
            gauge.appendChild(onlineAgentCount)
            free_workers = self.get_worker_count(connected=True, busy=False)
            busy_workers = self.get_worker_count(connected=True, busy=True)
            count = doc.createTextNode(str(free_workers + busy_workers))
            onlineAgentCount.appendChild(count)

            return doc, onlineAgentCount

        def add_offlineAgentCount(doc, gauge):
            offlineAgentCount = doc.createElement('offlineAgentCount')
            gauge.appendChild(offlineAgentCount)
            worker_count = self.get_worker_count(connected=False, busy=False)
            count = doc.createTextNode(str(worker_count))
            offlineAgentCount.appendChild(count)

            return doc, offlineAgentCount

        def add_workingAgentCount(doc, gauge):
            workingAgentCount = doc.createElement('workingAgentCount')
            gauge.appendChild(workingAgentCount)
            worker_count = self.get_worker_count(connected=True, busy=True)
            count = doc.createTextNode(str(worker_count))
            workingAgentCount.appendChild(count)

            return doc, workingAgentCount

        def add_availableAgentCount(doc, gauge):
            availableAgentCount = doc.createElement('availableAgentCount')
            gauge.appendChild(availableAgentCount)
            # worker_count = self.get_worker_count(connected=True)
            worker_count = self.get_worker_count(connected=True, busy=False)
            count = doc.createTextNode(str(worker_count))
            availableAgentCount.appendChild(count)

            return doc, availableAgentCount

        def add_unavailableAgentCount(doc, gauge):
            unavailableAgentCount = doc.createElement('unavailableAgentCount')
            gauge.appendChild(unavailableAgentCount)
            worker_count = self.get_worker_count(connected=True, busy=True)
            count = doc.createTextNode(str(worker_count))
            unavailableAgentCount.appendChild(count)

            return doc, unavailableAgentCount

        def add_workingMegaHertz(doc, gauge):
            workingMegaHertz = doc.createElement('workingMegaHertz')
            gauge.appendChild(workingMegaHertz)
            cpu_speed = self.get_cpu_speed(connected=True)
            mhz = doc.createTextNode(str(cpu_speed))
            workingMegaHertz.appendChild(mhz)

            return doc, workingMegaHertz

        def add_availableProcessorCount(doc, gauge):
            availableProcessorCount = doc.createElement(
                                                'availableProcessorCount')
            gauge.appendChild(availableProcessorCount)
            cpu_count = (self.get_cpu_count(connected=True) +
                         self.get_cpu_count(connected=False))

        def add_unavailableProcessorCount(doc, gauge):
            pass

        def add_onlineProcessorCount(doc, gauge):
            onlineProcessorCount = doc.createElement('onlineProcessorCount')
            gauge.appendChild(onlineProcessorCount)
            cpu_count = self.get_cpu_count(connected=True)
            c = doc.createTextNode(str(cpu_count))
            onlineProcessorCount.appendChild(c)

            return doc, onlineProcessorCount

        def add_offlineProcessorCount(doc, gauge):
            offlineProcessorCount = doc.createElement('offlineProcessorCount')
            gauge.appendChild(offlineProcessorCount)
            cpu_count = self.get_cpu_count(connected=False)
            c = doc.createTextNode(str(cpu_count))
            offlineProcessorCount.appendChild(c)

            return doc, offlineProcessorCount

        def add_workingProcessorCount(doc, gauge):
            workingProcessorCount = doc.createElement('workingProcessorCount')
            gauge.appendChild(workingProcessorCount)
            worker_count = self.get_cpu_count(connected=True)
            pcount = doc.createTextNode(str(worker_count))
            workingProcessorCount.appendChild(pcount)

            return doc, workingProcessorCount

        def add_workingAgentPercentage(doc, gauge):
            workingAgentPercentage = doc.createElement(
                                                    'workingAgentPercentage')
            gauge.appendChild(workingAgentPercentage)
            working_workers = self.get_worker_count(connected=True, busy=True)
            free_workers = self.get_worker_count(connected=True, busy=False)
            disconnected_workers = self.get_worker_count(connected=False,
                                                         busy=False)
            total_workers = (working_workers +
                             free_workers +
                             disconnected_workers)

            worker_percentage = float(working_workers / total_workers) * 100
            percentage = doc.createTextNode(str(worker_percentage))
            workingAgentPercentage.appendChild(percentage)

            return doc, workingAgentPercentage

        def add_date(doc, gauge):
            date = datetime.datetime.now()

            year = doc.createElement('Year')
            gauge.appendChild(year)
            year.appendChild(doc.createTextNode(str(date.year)))

            seconds = doc.createElement('Seconds')
            gauge.appendChild(seconds)
            seconds.appendChild(doc.createTextNode(str(date.second)))

            minutes = doc.createElement('Minutes')
            gauge.appendChild(minutes)
            minutes.appendChild(doc.createTextNode(str(date.minute)))

            return doc, year, seconds, minutes

        doc = xml.dom.minidom.Document()
        doc, gauge = create_gauge(doc)

        add_totalAgentCount(doc, gauge)
        add_onlineAgentCount(doc, gauge)
        add_offlineAgentCount(doc, gauge)
        add_availableAgentCount(doc, gauge)
        add_unavailableAgentCount(doc, gauge)
        add_workingAgentCount(doc, gauge)
        add_workingAgentPercentage(doc, gauge)

        add_onlineProcessorCount(doc, gauge)
        add_offlineAgentCount(doc, gauge)
        add_availableProcessorCount(doc, gauge)
        add_unavailableProcessorCount(doc, gauge)
        add_workingProcessorCount(doc, gauge)
        add_workingMegaHertz(doc, gauge)

        add_date(doc, gauge)
        s = cStringIO.StringIO()
        doc.writexml(s, newl='\n')

        return s.getvalue()