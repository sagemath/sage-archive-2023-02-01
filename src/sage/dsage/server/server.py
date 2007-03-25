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
##############################################################################

import zlib
import cPickle
import datetime
import xml
import cStringIO

from twisted.spread import pb
from twisted.python import log

from sage.dsage.misc.hostinfo import HostInfo
import sage.dsage.server.client_tracker as client_tracker
from sage.dsage.server.hostinfo_tracker import hostinfo_list
from sage.dsage.errors.exceptions import BadTypeError
from sage.dsage.database.job import expand_job

pb.setUnjellyableForClass(HostInfo, HostInfo)

class DSageServer(pb.Root):
    r"""
    This class represents Distributed Sage server which does all the
    coordination of distributing jobs, creating new jobs and accepting job
    submissions.

    """
    def __init__(self, jobdb, monitordb, log_level=0):
        r"""
        Initializes the Distributed Sage PB Server.

        Parameters:
        jobdb -- pass in the Database object
        log_level -- specifies the amount of logging to be done (default=0)

        """

        self.jobdb = jobdb
        self.monitordb = monitordb
        self.LOG_LEVEL = log_level

    def unpickle(self, pickled_job):
        return cPickle.loads(zlib.decompress(pickled_job))

    def get_job(self, anonymous=False, uuid=None):
        r"""
        Returns a job to the client.

        This method returns the first job that has not been completed
        in our job database.

        """

        if anonymous:
            jdict = self.jobdb.get_job(anonymous=True)
        else:
            jdict = self.jobdb.get_job(anonymous=False)
        if jdict == None:
            if self.LOG_LEVEL > 3:
                log.msg('[DSage, get_job]' + ' Job db is empty.')
            return None
        else:
            if self.LOG_LEVEL > 3:
                log.msg('[DSage, get_job]' + ' Returning Job %s to client'
                        % (jdict['job_id']))
            jdict['status'] = 'processing'
            self.jobdb.store_job(jdict)

        return jdict

    def set_job_uuid(self, job_id, uuid):
        return self.jobdb.set_job_uuid(job_id, uuid)

    def set_busy(self, uuid, busy):
        return self.monitordb.set_busy(uuid, busy=busy)

    def get_job_by_id(self, job_id):
        r"""
        Returns a job by the job id.

        Parameters:
        id -- the job id

        """

        job = self.jobdb.get_job_by_id(job_id)
        return job

    def get_job_result_by_id(self, job_id):
        """Returns the job result.

        Parameters:
        id -- the job id (str)

        """

        jdict = self.jobdb.get_job_by_id(job_id)
        job = expand_job(jdict)
        return job.result

    def get_job_output_by_id(self, job_id):
        """Returns the job output.

        Parameters:
        id -- the job id (str)

        """

        job = self.jobdb.get_job_by_id(job_id)

        return job.output

    def sync_job(self, job_id):
        job = self.jobdb.get_job_by_id(job_id)
        # new_job = copy.deepcopy(job)
        # print new_job
        # # Set file, data to 'Omitted' so we don't need to transfer it
        # new_job.code = 'Omitted...'
        # new_job.data = 'Omitted...'

        return job.pickle()

    def get_jobs_by_user_id(self, user_id):
        r"""
        Returns jobs created by user_id.

        Parameters:
        user_id -- the username (str)
        is_active -- when set to True, only return active jobs (bool)
        job_name -- the job name (optional)

        """

        jobs = self.jobdb.get_jobs_by_user_id(user_id)

        if self.LOG_LEVEL > 3:
            log.msg(jobs)
        return jobs

    def submit_job(self, jdict):
        r"""
        Submits a job to the job database.

        Parameters:
        jdict -- the internal dictionary of a Job object

        """


        if self.LOG_LEVEL > 3:
            log.msg('[DSage, submit_job] %s' % (jdict))

        if jdict['code'] is None:
            return False
        if jdict['name'] is None:
            jdict['name'] = 'No name specified'
        jdict['update_time'] = datetime.datetime.now()

        return self.jobdb.store_job(jdict)

    def get_all_jobs(self):
        r"""
        Returns a list of all jobs in the database.

        """
        return self.jobdb.get_all_jobs()

    def get_active_jobs(self):
        r"""
        Returns a list of active jobs"""

        return self.jobdb.get_active_jobs()

    def get_active_clients_list(self):
        r"""
        Returns a list of active clients.

        """

        raise NotImplementedError

    def get_killed_jobs_list(self):
        r"""
        Returns a list of killed job jdicts.
        """

        killed_jobs = self.jobdb.get_killed_jobs_list()
        return killed_jobs

    def get_next_job_id(self):
        r"""
        Returns the next job id.

        """

        if self.LOG_LEVEL > 0:
            log.msg('[DSage, get_next_job_id] Returning next job ID')

        return self.jobdb.get_next_job_id()

    def job_done(self, job_id, output, result, completed, worker_info):
        r"""
        job_done is called by the workers checkForJobOutput method.

        Parameters:
        job_id -- job id (str)
        output -- the stdout from the worker (string)
        result -- the result from the client (compressed pickle string)
                  result could be 'None'
        completed -- whether or not the job is completed (bool)
        worker_info -- ''.join(os.uname())

        """

        if self.LOG_LEVEL > 0:
            log.msg('[DSage, job_done] Job %s called back' % (job_id))
        if self.LOG_LEVEL > 3:
            log.msg('[DSage, job_done] Output: %s ' % output)
            log.msg('[DSage, job_done] Result: Some binary data...')
            log.msg('[DSage, job_done] completed: %s ' % completed)
            log.msg('[DSage, job_done] worker_info: %s ' % str(worker_info))

        jdict = self.get_job_by_id(job_id)

        if self.LOG_LEVEL > 3:
            log.msg('[DSage, job_done] result type' , type(result))

        output = str(output)
        if jdict['output'] is not None: # Append new output to existing output
            jdict['output'] += output
        else:
            jdict['output'] = output
        if completed:
            jdict['result'] = result
            jdict['status'] = 'completed'
            jdict['worker_info'] = str(worker_info)

        jdict['update_time'] = datetime.datetime.now()

        return self.jobdb.store_job(jdict)

    def job_failed(self, job_id):
        r"""
        job_failed is called when a remote job fails.

        Parameters:
        job_id -- the job id (str)

        """

        job = expand_job(self.jobdb.get_job_by_id(job_id))
        job.failures += 1

        if job.failures > self.jobdb.JOB_FAILURE_THRESHOLD:
            job.status = 'failed'
        else:
            job.status = 'new' # Put job back in the queue

        if self.LOG_LEVEL > 1:
            s = ['[DSage, job_failed] Job %s failed ' % (job_id),
                 '%s times. ' % (job.failures)]
            log.msg(''.join(s))

        return self.jobdb.store_job(job.reduce())

    def kill_job(self, job_id, reason):
        r"""
        Kills a job.

        Marks as job as killed and moves it to the killed jobs database.

        """

        if job_id == None:
            if self.LOG_LEVEL > 0:
                log.msg('[DSage, kill_job] No such job id %s' % job_id)
            return None
        else:
            self.jobdb.set_killed(job_id, killed=True)
            if self.LOG_LEVEL > 0:
                log.msg('Killed job %s' % (job_id))

        return job_id

    def get_monitor_list(self):
        r"""
        Returns a list of workers as a 3 tuple.

        tuple[0] = broker object
        tuple[1] = ip
        tuple[2] = port

        """
        return self.monitordb.get_monitor_list()
        # return worker_tracker.worker_list

    def get_client_list(self):
        r"""
        Returns a list of clients.

        """

        return client_tracker.client_list

    def get_cluster_speed(self):
        r"""
        Returns an approximation of the total CPU speed of the cluster.

        """

        cluster_speed = 0
        if self.LOG_LEVEL > 3:
            log.msg(hostinfo_list)
            log.msg(len(hostinfo_list))
        for h in hostinfo_list:
            speed_multiplier = int(h['cpus'])
            for k,v in h.iteritems():
                if k == 'cpu_speed':
                    cluster_speed += float(v) * speed_multiplier

        return cluster_speed

    def submit_host_info(self, h):
        r"""
        Takes a dict of workers machine specs.

        """

        if self.LOG_LEVEL > 0:
            log.msg(h)
        if len(hostinfo_list) == 0:
            hostinfo_list.append(h)
        else:
            for h in hostinfo_list:
                if h['uuid'] not in h.values():
                    hostinfo_list.append(h)

    def generate_xml_stats(self):
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
            working_workers = self.monitordb.get_worker_count(connected=True,
                                                              busy=True)
            free_workers = self.monitordb.get_worker_count(connected=True,
                                                           busy=False)
            disconnected_workers = self.monitordb.get_worker_count(
                                   connected=False,
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
            free_workers = self.monitordb.get_worker_count(connected=True,
                                                           busy=False)
            busy_workers = self.monitordb.get_worker_count(connected=True,
                                                           busy=True)
            count = doc.createTextNode(str(free_workers + busy_workers))
            onlineAgentCount.appendChild(count)

            return doc, onlineAgentCount

        def add_offlineAgentCount(doc, gauge):
            offlineAgentCount = doc.createElement('offlineAgentCount')
            gauge.appendChild(offlineAgentCount)
            worker_count = self.monitordb.get_worker_count(connected=False,
                                                           busy=False)
            count = doc.createTextNode(str(worker_count))
            offlineAgentCount.appendChild(count)

            return doc, offlineAgentCount

        def add_workingAgentCount(doc, gauge):
            workingAgentCount = doc.createElement('workingAgentCount')
            gauge.appendChild(workingAgentCount)
            worker_count = self.monitordb.get_worker_count(connected=True,
                                                           busy=True)
            count = doc.createTextNode(str(worker_count))
            workingAgentCount.appendChild(count)

            return doc, workingAgentCount

        def add_availableAgentCount(doc, gauge):
            availableAgentCount = doc.createElement('availableAgentCount')
            gauge.appendChild(availableAgentCount)
            # worker_count = self.monitordb.get_worker_count(connected=True)
            worker_count = self.monitordb.get_worker_count(connected=True,
                                                           busy=False)
            count = doc.createTextNode(str(worker_count))
            availableAgentCount.appendChild(count)

            return doc, availableAgentCount

        def add_unavailableAgentCount(doc, gauge):
            unavailableAgentCount = doc.createElement('unavailableAgentCount')
            gauge.appendChild(unavailableAgentCount)
            worker_count = self.monitordb.get_worker_count(connected=True,
                                                           busy=True)
            count = doc.createTextNode(str(worker_count))
            unavailableAgentCount.appendChild(count)

            return doc, unavailableAgentCount

        def add_workingMegaHertz(doc, gauge):
            workingMegaHertz = doc.createElement('workingMegaHertz')
            gauge.appendChild(workingMegaHertz)
            cpu_speed = self.monitordb.get_cpu_speed(connected=True)
            mhz = doc.createTextNode(str(cpu_speed))
            workingMegaHertz.appendChild(mhz)

            return doc, workingMegaHertz

        def add_availableProcessorCount(doc, gauge):
            pass

        def add_unavailableProcessorCount(doc, gauge):
            pass

        def add_onlineProcessorCount(doc, gauge):
            onlineProcessorCount = doc.createElement('onlineProcessorCount')
            gauge.appendChild(onlineProcessorCount)
            cpu_count = self.monitordb.get_cpu_count(connected=True)
            c = doc.createTextNode(str(cpu_count))
            onlineProcessorCount.appendChild(c)

            return doc, onlineProcessorCount

        def add_offlineProcessorCount(doc, gauge):
            offlineProcessorCount = doc.createElement('offlineProcessorCount')
            gauge.appendChild(offlineProcessorCount)
            cpu_count = self.monitordb.get_cpu_count(connected=False)
            c = doc.createTextNode(str(cpu_count))
            offlineProcessorCount.appendChild(c)

            return doc, offlineProcessorCount

        def add_workingProcessorCount(doc, gauge):
            workingProcessorCount = doc.createElement('workingProcessorCount')
            gauge.appendChild(workingProcessorCount)
            worker_count = self.monitordb.get_cpu_count(connected=True)
            pcount = doc.createTextNode(str(worker_count))
            workingProcessorCount.appendChild(pcount)

            return doc, workingProcessorCount

        def add_workingAgentPercentage(doc, gauge):
            workingAgentPercentage = doc.createElement(
                                                    'workingAgentPercentage')
            gauge.appendChild(workingAgentPercentage)
            working_workers = self.monitordb.get_worker_count(connected=True,
                                                              busy=True)
            free_workers = self.monitordb.get_worker_count(connected=True,
                                                           busy=False)
            disconnected_workers = self.monitordb.get_worker_count(
                                   connected=False,
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

class DSageWorkerServer(DSageServer):
    r"""
    Exposes methods to workers.
    """

    def remote_get_job(self):
        return DSageServer.get_job(self)

    def remote_job_done(self, job_id, output, result, completed, worker_info):
        if not (isinstance(job_id, str) or isinstance(completed, bool)):
            log.msg('BadType in remote_job_done')
            raise BadTypeError()

        return DSageServer.job_done(self, job_id, output, result,
                             completed, worker_info)

    def remote_job_failed(self, job_id):
        if not isinstance(job_id, str):
            log.msg('BadType in remote_job_failed')
            raise BadTypeError()

        return DSageServer.job_failed(self, job_id)

    def remote_get_killed_jobs_list(self):
        return DSageServer.get_killed_jobs_list(self)

    def remote_submit_host_info(self, hostinfo):
        if not isinstance(hostinfo, dict):
            log.msg('BadType in remote_submit_host_info')
            raise BadTypeError()
        return DSageServer.submit_host_info(self, hostinfo)

