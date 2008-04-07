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
import os
import sqlite3
from cStringIO import StringIO
import socket

from twisted.web2 import http, resource
from twisted.web2 import static
from twisted.conch.ssh.keys import Key
from twisted.conch.ssh.keys import BadKeyError

from sage.dsage.misc.constants import TMP_WORKER_FILES
from sage.dsage.server.stats import XMLStats

from xml.etree.ElementTree import (ElementTree as ET, Element, SubElement)

SAGE_ROOT  = os.environ['SAGE_ROOT']
DSAGE_LOCAL = SAGE_ROOT + '/local/dsage'
INDEX = os.path.join(DSAGE_LOCAL,'web/index.html')
STATIC = os.path.join(DSAGE_LOCAL,'web/static')

def create_jobs_table(jobs):
    """
    Returns an HTML table of jobs given a list of job dictionaries.

    """

    html = """
    """

    for i, job in enumerate(jobs):
        creation_time = job.creation_time.strftime('%F %r')
        try:
            update_time = job.update_time.strftime('%F %r')
        except AttributeError: # This is a fix for older databases.
            update_time = "N/A"
        html += """
        <tr class='tr%s'>
        """ % (i % 2)
        html += """
            <td><a href='#%s' onClick="getJobDetails('%s')">%s</a></td>
            <td class='job_status'>%s</td>
            <td>%s</td>
            <td>%s</td>
            <td>%s</td>
            <td>%s</td>
            <td>%s</td>
        </tr>
        """ % (job.job_id, job.job_id, job.name, job.status, job.username,
               creation_time, update_time, job.wall_time, job.priority)
    html += """"""

    return html

class Toplevel(resource.Resource):
    addSlash = True

    def __init__(self, dsage_server, web_port):
        self.dsage_server = dsage_server
        self.hostname = socket.getfqdn();
        self.index = open(INDEX).read() % (self.hostname,
                                           self.dsage_server.web_port,
                                           self.hostname,
                                           self.dsage_server.port,
                                           self.hostname,
                                           self.dsage_server.port)

    def child_static(self, ctx):
        return static.File(STATIC)

    def child_get_details(self, ctx):
        return GetJobDetails(self.dsage_server)

    def child_get_jobs(self, ctx):
        return GetJobs(self.dsage_server)

    def child_get_page(self, ctx):
        return GetPage(self.dsage_server)

    def child_pages(self, ctx):
        return Pages(self.dsage_server)

    def child_get_server_details(self, ctx):
        return GetServerDetails(self.dsage_server)

    def child_get_help(self, ctx):
        return GetHelp()

    def child_worker_files(self, ctx):
        return static.File(TMP_WORKER_FILES)

    def child_get_workers(self, ctx):
        return GetWorkers(self.dsage_server)

    def child_get_clients(self, ctx):
        return GetClients(self.dsage_server)

    def child_user_management(self, ctx):
        return UserManagement(self.dsage_server)

    def render(self, ctx):
        return http.Response(stream=self.index)

class GetHelp(resource.PostableResource):
    """
    Returns the help page.

    """

    def render(self, request):
        return static.File(os.path.join(STATIC, 'README.html'))

class Pages(resource.PostableResource):
    def __init__(self, dsage_server):
        self.dsage_server = dsage_server

    def render(self, request):
        count = int(request.args['count'][0])
        pages = self.dsage_server.jobdb.get_job_count() / count

        return http.Response(stream=str(pages))

class GetPage(resource.PostableResource):
    def __init__(self, dsage_server):
        self.dsage_server = dsage_server

    def render(self, request):
        page_num = int(request.args['n'][0])
        count = int(request.args['count'][0])
        job_count = self.dsage_server.jobdb.get_job_count()
        start = page_num * count;
        end = start + count;
        jobs = self.dsage_server.jobdb.get_job_range(start, end)

        html = create_jobs_table(jobs)

        return http.Response(stream=html)

class GetJobs(resource.PostableResource):
    """
    This resource returns a list of jobs in XML format to the client.

    """

    def __init__(self, dsage_server):
        self.dsage_server = dsage_server

    def render(self, request):
        try:
            count = int(request.args['count'][0])
        except:
            if count > 50:
                count = 50
            else:
                count = 10
        jobs = self.dsage_server.jobdb.get_n_jobs(count)

        html = create_jobs_table(jobs)

        return http.Response(stream=html)

    def build_xml(self, jdicts):
        """
        Builds an XML structure from a list of jdicts.

        """

        root = Element('jobs')

        for jdict in jdicts:
            job = SubElement(root, 'job')
            for k, v in jdict.iteritems():
                if k not in ('job_id', 'status', 'username', 'creation_time',
                             'last_update', 'priority'):
                    continue
                job.set(k, str(v))

        xml_stream = StringIO()
        tree = ET(root) # Wrap it into an ElementTree
        tree.write(xml_stream)

        return xml_stream

class GetJobDetails(resource.PostableResource):
    """
    This resource responds with details about a particular job.

    """

    def __init__(self, dsage_server):
        self.dsage_server = dsage_server

    def clean_code(self, code):
        if 'ans' in code:
            return code.split('\n')[0].split(' ')[2]
        else:
            return code

    def render(self, request):
        job_id = request.args['job_id'][0]
        try:
            html = """
            <thead>
            <tr>
                <th>Key</th>
                <th>Value</th>
            </tr>
            </thead>
            <tbody>
            """

            jdict = self.dsage_server.get_job_by_id(job_id)
            if not isinstance(jdict, dict):
                raise TypeError
            for i, (k, v) in enumerate(jdict.iteritems()):
                if k in ('data', 'id'):
                    continue
                if k == 'code': # We will display the code below the table
                    v = self.clean_code(v)
                if k == 'result': # result is an .sobj, link to the file
                    v = "<a href='worker_files/%s'>result.sobj</a>" % (jdict['job_id'] + '/' + 'result.sobj')
                html += """
                <tr class='tr%s'>
                """ % (i % 2)
                html += """
                <td>%s</td>
                <td>%s</td>
                </tr>
                """ % (k, v)

            html += """
            </tbody>
            """
        except (sqlite3.InterfaceError, TypeError):
            html = 'Invalid job id.'

        return http.Response(stream=html)

class GetWorkers(resource.PostableResource):
    """
    Returns a table of connected clients.

    """
    def __init__(self, dsage_server):
        self.dsage_server = dsage_server
        self.clientdb = self.dsage_server.clientdb
        self.workerdb = self.dsage_server.workerdb

    def gen_html(self):
        html = """
        <thead>
        <tr>
            <th colspan=9>Workers</th>
        </tr>
        <tr>
            <th>Username</th>
            <th>IP</th>
            <th># Workers</th>
            <th>Busy</th>
            <th>CPUs</th>
            <th>CPU Mhz</th>
            <th>Memory (MB)</th>
            <th>Connected</th>
            <th>Authenticated</th>
        </tr>
        <tbody>
        """

        for worker in self.workerdb.get_worker_list():
            html += """
            <tr>
                <td>%s</td>
                <td>%s</td>
                <td>%s</td>
                <td>%s</td>
                <td>%s</td>
                <td>%s</td>
                <td>%s</td>
                <td>%s</td>
                <td>%s</td>
            </tr>
            """ % (worker.username, worker.ip, worker.workers,worker.busy,
                   worker.cpus, worker.cpu_speed, worker.mem_total,
                   worker.connected, worker.authenticated)

        html += """
        </tbody>

        """

        return html

    def render(self, request):
        """
        Asks the server to generate stats and returns it in an XML file.

        """

        return http.Response(stream=self.gen_html())


class GetClients(resource.PostableResource):
    """
    Returns a table of connected clients.

    """
    def __init__(self, dsage_server):
        self.dsage_server = dsage_server
        self.clientdb = self.dsage_server.clientdb
        self.workerdb = self.dsage_server.workerdb

    def gen_html(self):
        html = """
        <thead>
        <tr>
            <th colspan=2>Clients</th>
        </tr>
        <tr>
            <th>Username</th>
            <th>Login time</th>
        </tr>
        <tbody>
        """

        for client in self.clientdb.get_client_list():
            html += """
            <tr>
                <td>%s</td>
                <td>%s</td>
            </tr>
            """ % (client.username, client.last_login.strftime('%F %T'))

        html += """
        </tbody>

        """

        return html

    def render(self, request):
        """
        Asks the server to generate stats and returns it in an XML file.

        """

        return http.Response(stream=self.gen_html())


class GetServerDetails(resource.PostableResource):
    """
    Returns an HTML table containing the server resources.

    """

    def __init__(self, dsage_server):
        self.dsage_server = dsage_server
        self.xml_stats = XMLStats(self.dsage_server)

    def gen_html(self):
        """
        generates html snippet from xml stats

        """

        self.xml_stats.gen_xml()

        html = """
        <thead>
        <tr>
        <th>Stat</th>
        <th>Value</th>
        </tr>
        <tbody>
        """

        for i, elem in enumerate(self.xml_stats.root.getchildren()):
            html += """
            <tr>
                <td>%s</td>
                <td>%s</td>
            </tr>
            """ % (' '.join(w.title() for w in elem.tag.split("_")),
                   elem.text)

        html += """
        </tbody>
        """

        return html

    def render(self, request):
        """
        Asks the server to generate stats and returns it in an XML file.

        """

        return http.Response(stream=self.gen_html())

class UserManagement(resource.PostableResource):
    """
    This page allows someone to add a user to the running server.

    """

    def __init__(self, dsage_server):
        self.dsage_server = dsage_server
        self.clientdb = dsage_server.clientdb

    def gen_html(self):
        html = """
        <form name="add_client" action="user_management" method="POST">
        Username:
        <input type="text" name="username">
        <br />
        Publickey:
        <br />
        <textarea id="pubkey" name="pubkey" cols=60 rows=10></textarea>
        <br />
        <input type="submit" value="Add User">
        </form>

        <form name="del_client" action="user_management" method="POST">
        Username:
        <input type=text" name="del_clientname">
        <br />
        <input type="submit" value="Delete User">
        <br />
        <br />
        """

        html += """
        <strong>Users</strong>
        """
        html += """<ul>"""
        for client in self.clientdb.get_client_list():
            html += """
            <li>%s</li>
            """ % client.username
        html += """
        </ul>
        """

        return html

    def render(self, request):
        if 'username' in request.args.keys():
            username = request.args['username'][0]
            try:
                pubkey = Key.fromString(request.args['pubkey'][0])
            except BadKeyError:
                return http.Response(stream="Invalid public key!")
            try:
                pubkey_str = pubkey.toString(type='openssh')
                self.clientdb.add_client(username, pubkey_str)
            except sqlite3.IntegrityError, msg:
                return http.Response(stream=msg)
        elif 'del_clientname' in request.args.keys():
            username = request.args['del_clientname'][0]
            try:
                self.clientdb.del_client(username)
            except Exception, msg:
                return http.Response(stream=msg)
            return http.Response(stream="User %s deleted!" % username)
        return http.Response(stream=self.gen_html())