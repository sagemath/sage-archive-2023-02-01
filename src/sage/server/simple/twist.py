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
        '{\n"a": [1,\n2,\n3],\n"b": "yep"\n}'
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
        return {
            'id': self.id,
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
