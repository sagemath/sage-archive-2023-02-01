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

from sage.dsage.misc.constants import SERVER_LOG
from sage.dsage.database.worker import Worker

class WorkerDatabaseSA(object):
    def __init__(self, Session):
        self.sess = Session()
        self._set_initial_state()

    def _set_initial_state(self):
        w = self.sess.query(Worker).all()
        for _w in w:
            _w.connected = False
            self.sess.save_or_update(_w)
        self.sess.commit()

    def set_authenticated(self, uuid, authenticated):
        w = self.sess.query(Worker).filter_by(uuid=uuid).first()
        w.authenticated = authenticated
        self.sess.save_or_update(w)
        self.sess.commit()

    def set_busy(self, uuid, busy):
        w = self.sess.query(Worker).filter_by(uuid=uuid).first()
        w.busy = busy
        self.sess.save_or_update(w)
        self.sess.commit()

    def add_worker(self, host_info):
        w = Worker(host_info)
        self.sess.save(w)
        self.sess.commit()

    def update_worker(self, host_info):
        uuid = host_info['uuid']
        w = self.sess.query(Worker).filter_by(uuid=uuid).first()
        for k, v in host_info.iteritems():
            setattr(w, k, v)
        self.sess.save_or_update(w)
        self.sess.commit()

    def get_worker(self, uuid):
        w = self.sess.query(Worker).filter_by(uuid=uuid).first()

        return w

    def get_worker_list(self):
        w = self.sess.query(Worker).all()

        return w

    def get_worker_by_job_id(self, job_id):
        w = self.sess.query(Worker).filter_by(job_id=job_id)

    def get_online_workers(self):
        w = self.sess.query(Worker).filter_by(connected=True).all()

        return w

    def get_worker_count(self, connected, busy):
        q = self.sess.query(Worker).filter_by(connected=connected,busy=busy)
        workers = q.all()

        count = sum([w.workers for w in workers])

        return count

    def get_cpu_speed(self, connected, busy):
        w = self.sess.query(Worker).all()

        return sum([_w.cpu_speed * _w.cpus for _w in w])

    def set_connected(self, uuid, connected):
        w = self.sess.query(Worker).filter_by(uuid=uuid).first()
        w.connected = connected
        self.sess.save_or_update(w)
        self.sess.commit()


class WorkerDatabase(object):
    """
    This table keeps track of workers.

    """

    def __init__(self, db_conn, log_file=SERVER_LOG, log_level=0):
        self.log_file = log_file
        self.log_level = log_level
        self.con = db_conn
        self.tablename = 'monitors'

    def _set_parameter(self, uuid, key, value):
        query = """UPDATE monitors
        SET %s=?
        WHERE uuid=?""" % (key)
        cur = self.con.cursor()
        cur.execute(query, (value, uuid))
        self.con.commit()

    def set_authenticated(self, uuid, authenticated):
        return self._set_parameter(uuid, 'authenticated', authenticated)

    def add_worker(self, host_info):
        query = """INSERT INTO monitors
        (uuid,
         username,
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
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """

        uuid = host_info['uuid']
        username = host_info['username']
        hostname = host_info['hostname']
        ip = host_info['ip']
        workers = host_info['workers']
        sage_version = host_info['sage_version']
        os_ = host_info['os']
        kernel_version = host_info['kernel_version']
        cpus = host_info['cpus']
        cpu_speed = host_info['cpu_speed']
        cpu_model = host_info['cpu_model']
        mem_total = host_info['mem_total']
        mem_free = host_info['mem_free']

        cur = self.con.cursor()
        cur.execute(query, (uuid, username, hostname, ip, workers,
                            sage_version, os_, kernel_version, cpus,
                            cpu_speed, cpu_model, mem_total, mem_free))
        self.con.commit()

    def update_worker(self, host_info):
        query = """UPDATE monitors
        SET hostname = ?, username = ?, ip = ?, workers = ?, sage_version = ?,
        os = ?, kernel_version = ?, cpus = ?, cpu_speed = ?, cpu_model = ?,
        mem_total = ?, mem_free = ? WHERE uuid = ?
        """

        uuid = host_info['uuid']
        username = host_info['username']
        hostname = host_info['hostname']
        ip = host_info['ip']
        workers = host_info['workers']
        sage_version = host_info['sage_version']
        os_ = host_info['os']
        kernel_version = host_info['kernel_version']
        cpus = host_info['cpus']
        cpu_speed = host_info['cpu_speed']
        cpu_model = host_info['cpu_model']
        mem_total = host_info['mem_total']
        mem_free = host_info['mem_free']

        cur = self.con.cursor()
        cur.execute(query, (hostname, username, ip, workers, sage_version,
                            os_, kernel_version, cpus, cpu_speed, cpu_model,
                            mem_total, mem_free, uuid))

    def get_worker(self, uuid):
        query = """SELECT
        uuid,
        workers,
        hostname,
        ip,
        authenticated,
        sage_version,
        os
        FROM monitors
        WHERE uuid = ?"""

        cur = self.con.cursor()
        cur.execute(query, (uuid,))
        result = cur.fetchone()
        if result is None:
            return result
        columns = [desc[0] for desc in cur.description]
        monitor = dict(zip(columns, result))
        for k, v in monitor.iteritems():
            if k == 'authenticated':
                monitor[k] = bool(v)

        return monitor

    def get_worker_list(self):
        """
        Returns a list of connected monitors.

        """

        query = """SELECT * FROM monitors"""
        cur = self.con.cursor()
        cur.execute(query)
        result = cur.fetchall()
        columns = [desc[0] for desc in cur.description]
        monitors = [dict(zip(columns, monitor)) for monitor in result]
        for monitor in monitors:
            for k, v in monitor.iteritems():
                 # Convert from 1/0 to python bool
                if k in ('authenticated', 'connected', 'busy'):
                    monitor[k] = bool(v)

        return monitors

    def set_connected(self, uuid, connected=True):
        """
        Sets the connected status of a monitor.

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

    def is_connected(self, uuid):
        """
        Returns whether the monitor is connected.

        """

        query = """SELECT connected FROM monitors WHERE uuid = ?"""
        cur = self.con.cursor()
        cur.execute(query, (uuid,))
        result = cur.fetchone()[0]

        return result

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

    def get_worker_count(self, connected, busy=False):
        """
        Returns the number of workers.

        Parameters:
        connected -- bool
        busy -- bool

        """

        if connected and not busy:
            query = """
            SELECT workers FROM monitors WHERE connected AND NOT busy
            """
        elif connected and busy:
            query = """
            SELECT workers FROM monitors WHERE connected AND busy
            """
        elif not connected and not busy:
            query = """
            SELECT workers FROM monitors WHERE NOT connected AND NOT busy
            """
        elif not connected and busy:
            query = """
            SELECT workers FROM monitors WHERE NOT connected AND busy
            """

        cur = self.con.cursor()
        cur.execute(query)

        result = cur.fetchall()

        return sum(w[0] for w in result)

    def get_cpu_speed(self, connected=True, busy=False):
        """
        Returns the aggregate cpu speed in Mhz.

        Parameters:
        connected -- bool

        """

        if connected and busy:
            query = """SELECT cpu_speed, workers FROM monitors
            WHERE connected AND busy"""
        elif connected:
            query = """SELECT cpu_speed, workers FROM monitors
            WHERE connected"""
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