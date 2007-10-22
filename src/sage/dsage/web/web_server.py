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

from twisted.web2 import http, resource
from twisted.web2 import static, http_headers, responsecode

SAGE_ROOT  = os.environ['SAGE_ROOT']
DSAGE_LOCAL = SAGE_ROOT + '/local/dsage'
CSS_FILE = os.path.join(DSAGE_LOCAL,'web/static/dsage_web.css')
SORTTABLE = os.path.join(DSAGE_LOCAL,'web/static/sorttable.js')
JS_FILE = os.path.join(DSAGE_LOCAL,'web/static/dsage_web.js')
PROTOTYPE = os.path.join(DSAGE_LOCAL,'web/static/prototype.js')
INDEX = os.path.join(DSAGE_LOCAL,'web/static/index.html')

# def create_jobs_table(jdicts):
#     """
#     Returns an HTML table of jobs given a list of job dictionaries.
#
#     """
#
#     html = """
#     <table class='jobs_table' cellspacing='0' padding='5px' border='solid #000 1px'>
#     <tr class='thead'>
#         <td>Job ID</td>
#         <td>Status</td>
#         <td>Username</td>
#         <td>Creation Time</td>
#         <td>Last Update</td>
#         <td>Priority</td>
#     </tr>
#     """
#
#     for i, jdict in enumerate(jdicts):
#         html+="""
#         <tr class='tr%s'
#         """ % (i % 2)
#         html+="""
#             <td><a href='get_details?job_id=%s'>%s</a></td>
#             <td>%s</td>
#             <td>%s</td>
#             <td>%s</td>
#             <td>%s</td>
#             <td>%s</td>
#         </tr>
#         """ % (jdict['job_id'],jdict['job_id'], jdict['status'],
#                jdict['username'], jdict['creation_time'],
#                jdict['update_time'], jdict['priority'])
#
#     html += '</table>'
#
#     return html

class Toplevel(resource.Resource):
    addSlash = True
    # header = """
    # <html>
    # <header>
    #     <link rel="stylesheet" type="text/css" href="dsage_web.css" />
    #     <scritp src='prototype.js' type='text/javascript'></script>
    #     <script src='dsage_web.js' type='text/javascript'></script>
    # </header>
    # <body>
    #     <div id='header'>
    #         <div id='title'><h1>DSAGE Job Status</h1></div>
    #     </div>
    # """
    #
    # footer = """
    # </body>
    # </html>
    # """

    def __init__(self, dsage_server):
        self.dsage_server = dsage_server
        self.child_get_details = GetJobDetails(self.dsage_server)
        self.child_get_jobs = GetJobs(self.dsage_server)

    def render(self, ctx):
        jobs = self.dsage_server.jobdb.get_all_jobs()
        # jobs_html = create_jobs_table(jobs)
        f = open(INDEX)
        return http.Response(stream=f.read())
        # return http.Response(stream=self.header + jobs_html + self.footer)
setattr(Toplevel, 'child_dsage_web.css', static.File(CSS_FILE))
setattr(Toplevel, 'child_dsage_web.js', static.File(JS_FILE))
setattr(Toplevel, 'child_prototype.js', static.File(PROTOTYPE))
setattr(Toplevel, 'child_sorttable.js', static.File(SORTTABLE))
# setattr(Toplevel, 'child_index.html', static.File(INDEX))

class GetJobs(resource.PostableResource):
    """
    This resource returns a list of jobs in XML format to the client.

    """

    def __init__(self, dsage_server):
        self.dsage_server = dsage_server
        self.jdicts = []

    def render(self, request):
        jdicts = self.dsage_server.jobdb.get_all_jobs()[:10]
        html = """
        <tr class='thead'>
            <td>Job ID</td>
            <td>Status</td>
            <td>Username</td>
            <td>Creation Time</td>
            <td>Last Update</td>
            <td>Priority</td>
        </tr>"""
        # if len(jdicts) != len(self.jdicts):
        #     self.jdicts = jdicts
        for i, jdict in enumerate(jdicts):
            html+="""
            <tr class='tr%s'>
            """ % (i % 2)
            html+="""
                <td><a href='get_details?job_id=%s'>%s</a></td>
                <td>%s</td>
                <td>%s</td>
                <td>%s</td>
                <td>%s</td>
                <td>%s</td>
            </tr>
            """ % (jdict['job_id'],jdict['job_id'], jdict['status'],
                   jdict['username'], jdict['creation_time'],
                   jdict['update_time'], jdict['priority'])

        html += '</table>'
        # else:
        #     pass
        # xml_stream = self.build_xml(jdicts)

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
        from cStringIO import StringIO
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

    header = """
    <html>
    <header>
        <link rel="stylesheet" type="text/css" href="dsage_web.css" />
        <scritp src='prototype.js' type='text/javascript'></script>
    </header>
    <body>
        <div id='header'>
            <div id='title'><h1>Details for %s</h1></div>
        </div>
    """

    footer = """
    </body>
    </html>
    """

    def __init__(self, dsage_server):
        self.dsage_server = dsage_server

    def render(self, request):
        job_id = request.args['job_id'][0]
        try:
            html = """
            <table class='job_details'>
            <tr class='thead'>
                <td>Key</td>
                <td>Value</td>
            </tr>
            """
            jdict = self.dsage_server.jobdb.get_job_by_id(job_id)
            if not isinstance(jdict, dict):
                raise TypeError
            for i, (k, v) in enumerate(jdict.iteritems()):
                if k == 'code': # We will display the code below the table
                    continue
                html += """
                <tr class='tr%s'>
                """ % (i % 2)
                html += """
                <td>%s</td>
                <td>%s</td>
                </tr>
                """ %(k, v)
            html += "</table>"
            html += """
            <div id='title'>
                Job Code
            </div>
            <div id='job_code'>
            %s
            <div>
            """ % jdict['code']
        except (sqlite3.InterfaceError, TypeError):
            html = 'Invalid job id.'
        html = self.header % (job_id) + html + self.footer
        return http.Response(stream=html)

class GetServerDetails(resource.Resource):
    """
    Returns an XML file containing the server resources.

    """

    def __init__(self, dsage_server):
        self.dsage_server = dsage_server

    def render(self):
        """
        Asks the server to generate stats and returns it in an XML file.

        """

        return self.dsage_server.generate_xml_stats()