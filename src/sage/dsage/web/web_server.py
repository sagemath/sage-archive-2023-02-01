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
from twisted.web2 import static, http_headers, responsecode

from sage.dsage.misc.constants import TMP_WORKER_FILES
from sage.dsage.server.stats import XMLStats

from xml.etree.ElementTree import (ElementTree as ET,
                                   Element,
                                   SubElement,
                                   dump,
                                   XML)

SAGE_ROOT  = os.environ['SAGE_ROOT']
DSAGE_LOCAL = SAGE_ROOT + '/local/dsage'
INDEX = os.path.join(DSAGE_LOCAL,'web/index.html')
STATIC = os.path.join(DSAGE_LOCAL,'web/static')

def create_jobs_table(jdicts):
    """
    Returns an HTML table of jobs given a list of job dictionaries.

    """

    html = """
    <thead>
    <tr>
        <th>Job ID</th>
        <th>Status</th>
        <th>Username</th>
        <th>Creation Time</th>
        <th>Last Update</th>
        <th>Priority</th>
    </tr>
    </thead>
    <tbody>
    """

    for i, jdict in enumerate(jdicts):
        creation_time = jdict['creation_time'].strftime('%F %r')
        try:
            update_time = jdict['update_time'].strftime('%F %r')
        except AttributeError: # This is a fix for older databases.
            update_time = "N/A"
        html+="""
        <tr class='tr%s'
        """ % (i % 2)
        html+="""
            <td><a href='#%s' onClick="getJobDetails('%s')">%s</a></td>
            <td>%s</td>
            <td>%s</td>
            <td>%s</td>
            <td>%s</td>
            <td>%s</td>
        </tr>
        """ % (jdict['job_id'], jdict['job_id'], jdict['job_id'],
               jdict['status'],
               jdict['username'], creation_time,
               update_time, jdict['priority'])
    html += """</tbody>"""

    return html

class Toplevel(resource.Resource):
    addSlash = True

    def __init__(self, dsage_server, server_port):
        self.dsage_server = dsage_server
        self.server_port = server_port

    def child_static(self, ctx):
        return static.File(STATIC)

    def child_get_details(self, ctx):
        return GetJobDetails(self.dsage_server)

    def child_get_jobs(self, ctx):
        return GetJobs(self.dsage_server)

    def child_get_server_details(self, ctx):
        return GetServerDetails(self.dsage_server)

    def child_get_help(self, ctx):
        return GetHelp()

    def child_worker_files(self, ctx):
        return static.File(TMP_WORKER_FILES)

    def render(self, ctx):
        index = open(INDEX).read() % (socket.getfqdn(), self.server_port)
        # return static.File(StringIO(index))
        return http.Response(stream=index)

class GetHelp(resource.PostableResource):
    """
    Returns the help page.

    """

    def render(self, request):
        return static.File(os.path.join(STATIC, 'README.html'))

class GetJobs(resource.PostableResource):
    """
    This resource returns a list of jobs in XML format to the client.

    """

    def __init__(self, dsage_server):
        self.dsage_server = dsage_server
        self.jdicts = []

    def render(self, request):
        try:
            count = int(request.args['count'][0])
        except:
            count = 10
        if count == - 1:
            jdicts = self.dsage_server.jobdb.get_all_jobs()
        else:
            jdicts = self.dsage_server.jobdb.get_all_jobs()[:count]

        html = create_jobs_table(jdicts)

        return http.Response(stream=html)

    def build_xml(self, jdicts):
        """
        Builds an XML structure from a list of jdicts.

        """


        from xml.etree.ElementTree import (ElementTree as ET,
                                           Element,
                                           SubElement,
                                           dump,
                                           XML)
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

            jdict = self.dsage_server.jobdb.get_job_by_id(job_id)
            if not isinstance(jdict, dict):
                raise TypeError
            for i, (k, v) in enumerate(jdict.iteritems()):
                # if k == 'code': # We will display the code below the table
                #     continue
                if k == 'result': # result is an .sobj, link to the file
                    v = "<a href='worker_files/%s'>result.sobj (Click to download.)</a>" % (jdict['job_id'] + '/' + 'result.sobj')
                html += """
                <tr class='tr%s'>
                """ % (i % 2)
                html += """
                <td>%s</td>
                <td>%s</td>
                </tr>
                """ %(k, v)

            html += """
            </tbody>
            """
        except (sqlite3.InterfaceError, TypeError):
            html = 'Invalid job id.'

        return http.Response(stream=html)

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
        <tbody>"""

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