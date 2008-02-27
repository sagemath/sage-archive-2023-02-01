r"""
This module provides a very simple API for interacting with a Sage session
over http. It runs in as part of the notebook server.

TESTS:
    sage: from sage.server.notebook.notebook_object import test_notebook
    sage: import urllib, re
    sage: def get_url(url): h = urllib.urlopen(url); data = h.read(); h.close(); return data

Start the notebook:
    sage: passwd = str(randint(1,1<<128))
    sage: nb = test_notebook(passwd, address='localhost', port=8095)

Login to a new session:
    sage: login_page = get_url('https://localhost:8095/simple/login?user=admin&password=%s' % passwd)
    sage: print login_page # random session info
    {
    "session": "2afee978c09b3d666c88b9b845c69608"
    }
    ___S_A_G_E___
    sage: session = re.match(r'.*"session": "([^"]*)"', login_page, re.DOTALL).groups()[0]

Run a command:
    sage: print get_url('https://localhost:8095/simple/compute?session=%s&code=2*2' % session)
    {
    "status": "done",
    "files": [],
    "cell_id": 1
    }
    ___S_A_G_E___
    4

Do a longer-running example:
    sage: n = next_prime(10^80)*next_prime(10^90)
    sage: print get_url('https://localhost:8095/simple/compute?session=%s&code=factor(%s)' % (session, n))
    {
    "status": "computing",
    "files": [],
    "cell_id": 2
    }
    ___S_A_G_E___

Get the status of the computation:
    sage: print get_url('https://localhost:8095/simple/status?session=%s&cell=2' % session)
    {
    "status": "computing",
    "files": [],
    "cell_id": 2
    }
    ___S_A_G_E___

Interrupt the computation:
    sage: _ = get_url('https://localhost:8095/simple/interrupt?session=%s' % session)

Test out getting files:
    sage: code = "h = open('a.txt', 'w'); h.write('test'); h.close()"
    sage: print get_url('https://localhost:8095/simple/compute?session=%s&code=%s' % (session, urllib.quote(code)))
    {
    "status": "done",
    "files": ["a.txt"],
    "cell_id": 3
    }
    ___S_A_G_E___

    sage: print get_url('https://localhost:8095/simple/file?session=%s&cell=3&file=a.txt' % session)
    test

Log out:
    sage: _ = get_url('https://localhost:8095/simple/logout?session=%s' % session)
    sage: nb.close(force=True)
"""

#############################################################################
#   Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################


import re, random, os.path

from twisted.internet.threads import deferToThread
from twisted.python import log
from twisted.internet import defer, reactor
from twisted.web.server import NOT_DONE_YET

from twisted.web2 import server, http, resource, channel
from twisted.web2 import static, http_headers, responsecode

sessions = {}

# There must be a better way to avoid circular imports...
late_import_done = False

def late_import():
    global SEP, notebook_twist, late_import_done
    if not late_import_done:
        from sage.server.notebook.twist import SEP
        import sage.server.notebook.twist as notebook_twist
        late_import_done = True

def simple_jsonize(data):
    """
    This may be replaced by a JSON spkg."

    EXAMPLES:
        sage: from sage.server.simple.twist import simple_jsonize
        sage: print simple_jsonize({'a': [1,2,3], 'b': "yep"})
        { "a": [1, 2, 3], "b": "yep" }
    """
    if isinstance(data, dict):
        values = ['"%s": %s' % (key, simple_jsonize(value)) for key, value in data.iteritems()]
        return "{\n%s\n}" % ',\n'.join(values)
    elif isinstance(data, (list, tuple)):
        values = [simple_jsonize(value) for value in data]
        return "[%s]" % ",\n".join(values)
    elif isinstance(data, bool):
        return str(data).lower()
    elif data is None:
        return "null"
    else:
        value =  str(data)
        if re.match(r'^\d+(\.\d*)?$', value):
            return value
        else:
            return '"%s"' % value

class SessionObject:
    def __init__(self, id, username, worksheet, timeout=1):
        self.id = id
        self.username = username
        self.worksheet = worksheet
        self.default_timeout = timeout

    def get_status(self):
        """
        Return a dictionary to be returned (in JSON format) representing
        the status of self.

        TEST:
            sage: from sage.server.simple.twist import SessionObject
            sage: s = SessionObject(id=1, username=None, worksheet=None)
            sage: s.get_status()
            {'session': 1}
        """
        return {
            'session': self.id
        }

class LoginResource(resource.Resource):

    def render(self, ctx):
        late_import()
        username = "admin" # ctx.args['username'][0]
        worksheet = notebook_twist.notebook.create_new_worksheet("_simple_session", username)
        worksheet.sage() # create the sage session
        worksheet.initialize_sage()
        # is this a secure enough random number?
        session_id = "%x" % random.randint(1, 1 << 128)
        session = SessionObject(session_id, username, worksheet)
        sessions[session_id] = session
        # FOR TESTING
        sessions['test'] = session
        status = session.get_status()
        return http.Response(stream = "%s\n%s\n" % (simple_jsonize(status), SEP))

class LogoutResource(resource.Resource):

    def render(self, ctx):
        late_import()
        session = sessions[ctx.args['session'][0]]
        session.worksheet.notebook().delete_worksheet(session.worksheet.filename())
        status = session.get_status()
        return http.Response(stream = "%s\n%s\n" % (simple_jsonize(status), SEP))

class InterruptResource(resource.Resource):

    def render(self, ctx):
        late_import()
        session = sessions[ctx.args['session'][0]]
        session.worksheet.interrupt()
        status = session.get_status()
        return http.Response(stream = "%s\n%s\n" % (simple_jsonize(status), SEP))

class RestartResource(resource.Resource):

    def render(self, ctx):
        late_import()
        session = sessions[ctx.args['session'][0]]
        session.worksheet.restart_sage()
        status = session.get_status()
        return http.Response(stream = "%s\n%s\n" % (simple_jsonize(status), SEP))

class SessionStatus:

    def render_session_status(self, session):
        status = {}
        return http.Response(stream = "\n".join([simple_jsonize(status), SEP]))

class CellResource(resource.Resource):

    def render_cell_result(self, cell):

        late_import()
        cell.worksheet().check_comp()
        if cell.interrupted():
            cell_status = 'interrupted'
        elif cell.computing():
            cell_status = 'computing'
        else:
            cell_status = 'done'
        status = { 'cell_id': cell.id(), 'status': cell_status, 'files': cell.files() }
        result = cell.output_text(raw=True)
        return http.Response(stream = "\n".join([simple_jsonize(status), SEP, result]))


class StatusResource(CellResource):

    def render(self, ctx):
        session = sessions[ctx.args['session'][0]]
        try:
            cell_id = int(ctx.args['cell'][0])
            cell = session.worksheet.get_cell_with_id(cell_id)
            return self.render_cell_result(cell)
        except KeyError:
            status = session.get_status()
            return http.Response(stream = "%s\n%s\n" % (simple_jsonize(status), SEP))

class FileResource(resource.Resource):
    """
    This differs from the rest as it does not print a header, just the raw file data.
    """
    def render(self, ctx):
        session = sessions[ctx.args['session'][0]]
        cell_id = int(ctx.args['cell'][0])
        cell = session.worksheet.get_cell_with_id(cell_id)
        file_name = ctx.args['file'][0]
        if file_name in cell.files():
            return static.File("%s/%s" % (cell.directory(), file_name))
        else:
            return http.Response(code=404, stream = "No such file %s in cell %s." % (file_name, cell_id))

class ComputeResource(CellResource):

    def render(self, ctx):
        late_import()
        session = sessions[ctx.args['session'][0]]
        # is this a secure enough random number?
        cell = session.worksheet.append_new_cell()
        cell.set_input_text(ctx.args['code'][0])
        cell.evaluate(username = session.username)
        # session.worksheet.start_next_comp()

        d = defer.Deferred()
        d.addCallback(self.render_cell_result)
        reactor.callLater(session.default_timeout, d.callback, cell) # is there a way to call this earlier?

        return d


class SimpleServer(resource.Resource):

    child_login = LoginResource()
    child_logout = LogoutResource()
    child_interrupt = InterruptResource()
    child_restart = RestartResource()

    child_compute = ComputeResource()
    child_status = StatusResource()
    child_file = FileResource()

    def render(self, ctx):
        return http.Response(stream="Yo!")
