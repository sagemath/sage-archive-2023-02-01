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
JS_FILE = os.path.join(DSAGE_LOCAL,'web/static/dsage_web.js')
PROTOTYPE = os.path.join(DSAGE_LOCAL,'web/static/prototype.js')
INDEX = os.path.join(DSAGE_LOCAL,'web/static/index.html')

def create_jobs_table(jdicts):
    """
    Returns an HTML table of jobs given a list of job dictionaries.

    """

    html = """
    <table class='jobs_table' cellspacing='0' padding='5px' border='solid #000 1px'>
    <tr class='thead'>
        <td>Job ID</td>
        <td>Status</td>
        <td>Username</td>
        <td>Creation Time</td>
        <td>Last Update</td>
        <td>Priority</td>
    </tr>
    """

    for jdict in jdicts:
        html+="""
        <tr class='tr0'>
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

    return html

class Toplevel(resource.Resource):
    addSlash = True
    header = """
    <html>
    <header>
        <link rel="stylesheet" type="text/css" href="dsage_web.css" />
        <scritp src='prototype.js' type='text/javascript'></script>
        <script src='dsage_web.js' type='text/javascript'></script>
    </header>
    <body>
        <div id='header'>
            <div id='title'><h1>DSAGE Job Status</h1></div>
        </div>
    """

    footer = """
    </body>
    </html>
    """

    def __init__(self, dsage_server):
        self.dsage_server = dsage_server
        self.child_get_details = GetJobDetails(self.dsage_server)

    def render(self, ctx):
        jobs = self.dsage_server.jobdb.get_all_jobs()
        jobs_html = create_jobs_table(jobs)
        f = open(INDEX)
        return http.Response(stream=f.read())
        # return http.Response(stream=self.header + jobs_html + self.footer)
setattr(Toplevel, 'child_dsage_web.css', static.File(CSS_FILE))
setattr(Toplevel, 'child_dsage_web.js', static.File(JS_FILE))
setattr(Toplevel, 'child_prototype.js', static.File(PROTOTYPE))
# setattr(Toplevel, 'child_index.html', static.File(INDEX))

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
            for k, v in jdict.iteritems():
                if k == 'code': # We will display the code below the table
                    continue
                html += """
                <tr class='tr0'>
                <td>%s</td>
                <td>%s</td>
                </tr>
                """ %(k, v)
            html += "</table>"
            html += """
            <div id='codeh'>
                Job Code
            </div>
            <div id='code'>
            %s
            <div>
            """ % jdict['code']
        except (sqlite3.InterfaceError, TypeError):
            html = 'Invalid job id.'
        return http.Response(stream=self.header % (job_id) + html + self.footer)