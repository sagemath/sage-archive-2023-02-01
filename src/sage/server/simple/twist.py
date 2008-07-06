r"""
This module provides a very simple API for interacting with a Sage session
over http. It runs in as part of the notebook server.

NOTE:
    The exact data in the JSON header may vary over time (for example, further
    data may be added), but should remain backwards compatable if it is being
    parsed as JSON data.


TESTS:
Start the notebook:
    sage: from sage.server.misc import find_next_available_port
    sage: port = find_next_available_port(9000, verbose=False)
    sage: from sage.server.notebook.notebook_object import test_notebook
    sage: passwd = str(randint(1,1<<128))
    sage: nb = test_notebook(passwd, secure=False, address='localhost', port=port, verbose=True) #doctest: +ELLIPSIS
    ...
    Notebook started.

Import urllib:
    sage: import urllib, re
    sage: def get_url(url): h = urllib.urlopen(url); data = h.read(); h.close(); return data


Login to a new session:
    sage: login_page = get_url('http://localhost:%s/simple/login?username=admin&password=%s' % (port, passwd))
    sage: print login_page # random session id
    {
    "session": "2afee978c09b3d666c88b9b845c69608"
    }
    ___S_A_G_E___
    sage: session = re.match(r'.*"session": "([^"]*)"', login_page, re.DOTALL).groups()[0]

Run a command:
    sage: sleep(0.5)
    sage: print get_url('http://localhost:%s/simple/compute?session=%s&code=2*2' % (port, session))
    {
    "status": "done",
    "files": [],
    "cell_id": 1
    }
    ___S_A_G_E___
    4

Do a longer-running example:
    sage: n = next_prime(10^25)*next_prime(10^30)
    sage: print get_url('http://localhost:%s/simple/compute?session=%s&code=factor(%s)&timeout=0.1' % (port, session, n))
    {
    "status": "computing",
    "files": [],
    "cell_id": 2
    }
    ___S_A_G_E___

Get the status of the computation:
    sage: print get_url('http://localhost:%s/simple/status?session=%s&cell=2' % (port, session))
    {
    "status": "computing",
    "files": [],
    "cell_id": 2
    }
    ___S_A_G_E___

Interrupt the computation:
    sage: _ = get_url('http://localhost:%s/simple/interrupt?session=%s' % (port, session))

Test out getting files:
    sage: code = "h = open('a.txt', 'w'); h.write('test'); h.close()"
    sage: print get_url('http://localhost:%s/simple/compute?session=%s&code=%s' % (port, session, urllib.quote(code)))
    {
    "status": "done",
    "files": ["a.txt"],
    "cell_id": 3
    }
    ___S_A_G_E___

    sage: print get_url('http://localhost:%s/simple/file?session=%s&cell=3&file=a.txt' % (port, session))
    test

Log out:
    sage: _ = get_url('http://localhost:%s/simple/logout?session=%s' % (port, session))
    sage: nb.dispose()
"""

#############################################################################
#   Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################


import re, random, os.path, shutil, time

from twisted.internet.task import LoopingCall
from twisted.python import log
from twisted.internet import defer, reactor
from twisted.cred import credentials

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
    def __init__(self, id, username, worksheet, timeout=5):
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
        try:
            username = ctx.args['username'][0]
            password = ctx.args['password'][0]
            U = notebook_twist.notebook.user(username)
            if not U.password_is(password):
                raise ValueError # want to return the same error message
        except KeyError, ValueError:
            return http.Response(stream = "Invalid login.")

        worksheet = notebook_twist.notebook.create_new_worksheet("_simple_session_", username, add_to_list=False)
        worksheet.sage() # create the sage session
        worksheet.initialize_sage()
        # is this a secure enough random number?
        session_id = "%x" % random.randint(1, 1 << 128)
        session = SessionObject(session_id, username, worksheet)
        sessions[session_id] = session
        sessions['test'] = session
        status = session.get_status()
        return http.Response(stream = "%s\n%s\n" % (simple_jsonize(status), SEP))

class LogoutResource(resource.Resource):

    def render(self, ctx):
        try:
            session = sessions[ctx.args['session'][0]]
        except KeyError:
            return http.Response(stream = "Invalid session.")
        session.worksheet.quit()
        shutil.rmtree(session.worksheet.directory())
        status = session.get_status()
        return http.Response(stream = "%s\n%s\n" % (simple_jsonize(status), SEP))

class InterruptResource(resource.Resource):

    def render(self, ctx):
        try:
            session = sessions[ctx.args['session'][0]]
        except KeyError:
            return http.Response(stream = "Invalid session.")
        session.worksheet.interrupt()
        status = session.get_status()
        return http.Response(stream = "%s\n%s\n" % (simple_jsonize(status), SEP))

class RestartResource(resource.Resource):

    def render(self, ctx):
        try:
            session = sessions[ctx.args['session'][0]]
        except KeyError:
            return http.Response(stream = "Invalid session.")
        session.worksheet.restart_sage()
        status = session.get_status()
        return http.Response(stream = "%s\n%s\n" % (simple_jsonize(status), SEP))

class CellResource(resource.Resource):

    def start_comp(self, cell, timeout):
        start_time = time.time()
        looper_list = []
        looper = LoopingCall(self.check_comp, cell, start_time, timeout, looper_list)
        looper_list.append(looper) # so check_comp has access
        looper.cell = cell # to pass it on
        d = looper.start(0.25, now=True)
        d.addCallback(self.render_cell_result)
        return d

    def check_comp(self, cell, start_time, timeout, looper_list):
        cell.worksheet().check_comp(wait=0.01) # don't want to block, delay handled by twisted
        if not cell.computing() or time.time() - start_time > timeout:
            looper_list[0].stop()

    def render_cell_result(self, looper):
        cell = looper.cell
        if cell.interrupted():
            cell_status = 'interrupted'
        elif cell.computing():
            cell_status = 'computing'
        else:
            cell_status = 'done'
        status = { 'cell_id': cell.id(), 'status': cell_status, 'files': cell.files() }
        result = cell.output_text(raw=True)
        return http.Response(stream = "\n".join([simple_jsonize(status), SEP, result]))


class ComputeResource(CellResource):

    def render(self, ctx):
        try:
            session = sessions[ctx.args['session'][0]]
        except KeyError:
            return http.Response(stream = "Invalid session.")
        try:
            timeout = float(ctx.args['timeout'][0])
        except (KeyError, ValueError), msg:
            timeout = session.default_timeout
        cell = session.worksheet.append_new_cell()
        cell.set_input_text(ctx.args['code'][0])
        cell.evaluate(username = session.username)
        return self.start_comp(cell, timeout)


class StatusResource(CellResource):

    def render(self, ctx):
        try:
            session = sessions[ctx.args['session'][0]]
        except KeyError:
            return http.Response(stream = "Invalid session.")
        try:
            cell_id = int(ctx.args['cell'][0])
            cell = session.worksheet.get_cell_with_id(cell_id)
            try:
                timeout = float(ctx.args['timeout'][0])
            except KeyError, ValueError:
                timeout = -1
            return self.start_comp(cell, timeout)
        except KeyError:
            status = session.get_status()
            return http.Response(stream = "%s\n%s\n" % (simple_jsonize(status), SEP))


class FileResource(resource.Resource):
    """
    This differs from the rest as it does not print a header, just the raw file data.
    """
    def render(self, ctx):
        try:
            session = sessions[ctx.args['session'][0]]
        except KeyError:
            return http.Response(code=404, stream = "Invalid session.")
        try:
            cell_id = int(ctx.args['cell'][0])
            file_name = ctx.args['file'][0]
        except KeyError:
            return http.Response(code=404, stream = "Unspecified cell or file.")
        cell = session.worksheet.get_cell_with_id(cell_id)
        if file_name in cell.files():
            return static.File("%s/%s" % (cell.directory(), file_name))
        else:
            return http.Response(code=404, stream = "No such file %s in cell %s." % (file_name, cell_id))


class SimpleServer(resource.Resource):

    child_login = LoginResource()
    child_logout = LogoutResource()
    child_interrupt = InterruptResource()
    child_restart = RestartResource()

    child_compute = ComputeResource()
    child_status = StatusResource()
    child_file = FileResource()

    def render(self, ctx):
        return http.Response(stream="Please login.")
