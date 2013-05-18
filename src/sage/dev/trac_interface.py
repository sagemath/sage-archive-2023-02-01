"""
Trac Interface
"""
import os
import re
import tempfile
import time
import urllib
import urllib2
import urlparse

from xmlrpclib import SafeTransport, ServerProxy

from sage.doctest import DOCTEST_MODE
from sage.env import REALM, TRAC_SERVER_URI

FIELD_REGEX = re.compile("^([A-Za-z ]+):(.*)$")
ALLOWED_FIELDS = {
        "authors":          "Authors",
        "branch":           "Branch",
        "cc":               "Cc",
        "component":        "Component",
        "dependencies":     "Dependencies",
        "keywords":         "Keywords",
        "merged in":        "Merged in",
        "milestone":        "Milestone",
        "owned by":         "Owned by",
        "priority":         "Priority",
        "report upstream":  "Report Upstream",
        "reviewers":        "Reviewers",
        "stopgaps":         "Stopgaps",
        "type":             "Type",
        "work issues":      "Work issues",
        }
TICKET_FILE_GUIDE = """
# Lines starting with `#` are ignored.
# Lines starting with `Field: ` correspond to fields of
# the trac ticket, and can be followed by text on the same line.
# They will be assigned to the corresponding field on the trac
# ticket.
#
# Lines not following this format will be put into the ticket
# description. Trac markup is supported.
#
# An empty file aborts ticket creation/editing.
"""

class TicketSyntaxError(SyntaxError): # we don't want to catch normal syntax errors
    pass

def _parse_ticket_file(filename):
    r"""
    parses ticket file

    INPUT:

    - ``filename`` -- filename to parse

    OUTPUT:

    - ``summary`` -- a string consisting of the ticket's summary

    - ``description`` -- a string consisting of the ticket's description

    - ``fields`` -- a dictionary containing field, entry pairs

    TESTS::

        sage: from sage.dev.trac_interface import _parse_ticket_file
        sage: import os, tempfile
        sage: tmp = tempfile.mkstemp()[1]
        sage: with open(tmp, 'w') as f:
        ....:     f.write("no summary\n")
        sage: _parse_ticket_file(tmp)
        Traceback (most recent call last):
        ...
        TicketSyntaxError: no valid summary found
        sage: with open(tmp, 'w') as f:
        ....:     f.write("summary:no description\n")
        sage: _parse_ticket_file(tmp)
        Traceback (most recent call last):
        ...
        TicketSyntaxError: no description found
        sage: with open(tmp, 'w') as f:
        ....:     f.write("summary:double summary\n")
        ....:     f.write("summary:double summary\n")
        sage: _parse_ticket_file(tmp)
        Traceback (most recent call last):
        ...
        TicketSyntaxError: line 2: only one value for summary allowed
        sage: with open(tmp, 'w') as f:
        ....:     f.write("bad field:bad field entry\n")
        sage: _parse_ticket_file(tmp)
        Traceback (most recent call last):
        ...
        TicketSyntaxError: line 1: field `bad field` not supported
        sage: with open(tmp, 'w') as f:
        ....:     f.write("summary:a summary\n")
        ....:     f.write("branch:a branch\n")
        ....:     f.write("some description\n")
        ....:     f.write("#an ignored line\n")
        ....:     f.write("more description\n")
        ....:     f.write("\n")
        sage: _parse_ticket_file(tmp)
        ('a summary', 'some description\nmore description\n', {'branch': 'a branch'})
        sage: with open(tmp, 'w') as f:
        ....:     f.write("summary:a summary\n")
        ....:     f.write("some description\n")
        ....:     f.write("branch:a branch\n")
        sage: _parse_ticket_file(tmp)
        ('a summary', 'some description\nbranch:a branch\n', {})
        sage: os.unlink(tmp)
    """
    lines = list(open(filename).read().splitlines())

    if all(l.rstrip().startswith('#') for l in lines if l.rstrip()):
        return

    fields = {}
    for i, line in enumerate(lines):
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        m = FIELD_REGEX.match(line)
        if m:
            field = m.groups()[0]
            if not (field.lower() == 'summary' or
                    field.lower() in ALLOWED_FIELDS):
                raise TicketSyntaxError("line %s: "%(i+1) +
                                        "field `%s` not supported"%field)
            elif field.lower() in fields:
                raise TicketSyntaxError("line %s: "%(i+1) +
                                        "only one value for %s allowed"%field)
            else:
                fields[field.lower()] = m.groups()[1].strip()
                continue
        else:
            break
    else: # no description
        i += 1

    # separate summary from other fields
    try:
        summary = fields.pop('summary')
    except KeyError:
        summary = None

    description = [line.rstrip() for line in lines[i:]
            if not line.startswith('#')]

    # remove leading and trailing empty newlines
    while description and not description[0]:
        description.pop(0)
    while description and not description[-1]:
        description.pop()

    if not summary:
        raise TicketSyntaxError("no valid summary found")
    elif not description:
        raise TicketSyntaxError("no description found")
    else:
        return summary, "\n".join(description)+"\n", fields

class DigestTransport(object, SafeTransport):
    """
    Handles an HTTP transaction to an XML-RPC server.

    EXAMPLES::

        sage: sage.dev.trac_interface.DigestTransport()
        <sage.dev.trac_interface.DigestTransport object at ...>
    """
    def __init__(self, **kwds):
        """
        Initialization.

        EXAMPLES::

            sage: type(sage.dev.trac_interface.DigestTransport())
            <class 'sage.dev.trac_interface.DigestTransport'>
            sage: type(sage.dev.trac_interface.DigestTransport(realm='realm',
            ....:         url='url', username='username', password='password'))
            <class 'sage.dev.trac_interface.DigestTransport'>
        """
        def get_pop(this, k, d=None):
            try:
                return this.pop(k)
            except KeyError:
                return d

        auth = tuple(get_pop(kwds, x) for x in
                ('realm', 'url', 'username', 'password'))

        SafeTransport.__init__(self, **kwds)

        authhandler = urllib2.HTTPDigestAuthHandler()
        if all(x is not None for x in auth):
            authhandler.add_password(*auth)

        self.opener = urllib2.build_opener(authhandler)

    def single_request(self, host, handler, request_body, verbose):
        """
        Issue an XML-RPC request.

        EXAMPLES::

            sage: from sage.env import TRAC_SERVER_URI
            sage: import urlparse
            sage: url = urlparse.urlparse(TRAC_SERVER_URI).netloc
            sage: d = sage.dev.trac_interface.DigestTransport()
            sage: d.single_request(url, 'xmlrpc',
            ....: '''<?xml version='1.0'?>
            ....: <methodCall>
            ....: <methodName>ticket.get</methodName>
            ....: <params>
            ....: <param>
            ....: <value><int>1000</int></value>
            ....: </param>
            ....: </params>
            ....: </methodCall>
            ....: ''', 0) # optional: online
            ([1000,
              <DateTime '20071025T16:48:05' at ...>,
              <DateTime '20080110T08:28:40' at ...>,
              {'status': 'closed',
               'changetime': <DateTime '20080110T08:28:40' at ...>,
               'description': '',
               'reporter': 'was',
               'cc': '',
               'type': 'defect',
               'milestone': 'sage-2.10',
               '_ts': '1199953720000000',
               'component': 'distribution',
               'summary': 'Sage does not have 10000 users yet.',
               'priority': 'major',
               'owner': 'was',
               'time': <DateTime '20071025T16:48:05' at ...>,
               'keywords': '',
               'resolution': 'fixed'}],)
        """
        req = urllib2.Request(
                urlparse.urlunparse(('http', host, handler, '', '', '')),
                request_body, {'Content-Type': 'text/xml',
                    'User-Agent': self.user_agent})

        response = self.opener.open(req)

        self.verbose = verbose
        return self.parse_response(response)

class DoctestServerProxy(object):
    """
    A fake trac proxy for doctesting the functionality in this file which would require authentication by trac.

    EXAMPLES::

        sage: sage.dev.trac_interface.DoctestServerProxy(None)
        <sage.dev.trac_interface.DoctestServerProxy object at ...>
    """
    def __init__(self, trac):
        self._trac = trac
        self._sshkeys = {}

    @property
    def sshkeys(self):
        try:
            return self._sshkeys_impl
        except AttributeError:
            pass
        class SshKeys(object):
            def setkeys(this, keys):
                if self._trac._username not in self._sshkeys:
                    self._sshkeys[self._trac._username] = set()
                self._sshkeys[self._trac._username] = self._sshkeys[self._trac._username].union(set(keys))
                return 0
            def getkeys(this):
                if self._trac._username not in self._sshkeys: return []
                return list(self._sshkeys[self._trac._username])
            def listusers(this):
                return self._sshkeys.keys()

        self._sshkeys_impl = SshKeys()
        return self._sshkeys_impl

    @property
    def ticket(self):
        class Ticket(object):
            def create(self, summary, description, attributes, notify):
                return 14366
            def update(self, ticketnum, comment, attributes):
                assert isinstance(ticketnum, int)
                return (ticketnum,)

        return Ticket()

class TracInterface(object):
    """
    Wrapper around the XML-RPC interface of trac.

    EXAMPLES::

        sage: dev.trac
        <sage.dev.trac_interface.TracInterface object at ...>
    """
    def __init__(self, sagedev):
        """
        Initialization.

        EXAMPLES::

            sage: type(sage.dev.trac_interface.TracInterface(dev))
            <class 'sage.dev.trac_interface.TracInterface'>
        """
        self._sagedev = sagedev
        self._UI = sagedev._UI
        sagedev._config.setdefault('trac', {})
        self._config = sagedev._config['trac']

    @property
    def _anonymous_server_proxy(self):
        """
        Lazy wrapper around a non-authenticated XML-RPC interface to trac.

        .. NOTE::

            Unlike the authenticated server proxy, this is not replaced with a
            fake proxy for doctesting. All doctests using it should therefore
            be labeled as optional ``online``

        EXAMPLES::

            sage: from sage.dev.sagedev import SageDev, doctest_config
            sage: conf  = doctest_config()
            sage: conf['trac']['server'] = 'http://trac.sagemath.org/'
            sage: trac = SageDev(conf).trac
            sage: trac._anonymous_server_proxy
            <ServerProxy for trac.sagemath.org/xmlrpc>
        """
        try:
            return self.__anonymous_server_proxy
        except AttributeError:
            pass

        realm = self._config.get('realm', REALM)
        server = self._config.get('server', TRAC_SERVER_URI)

        url = urlparse.urljoin(server, 'xmlrpc')

        transport = DigestTransport()
        self.__anonymous_server_proxy = ServerProxy(url, transport=transport)
        return self.__anonymous_server_proxy

    @property
    def _authenticated_server_proxy(self):
        """
        Get an XML-RPC proxy object that is authenticated using the users
        username and password.

        EXAMPLES::

            sage: dev.trac._authenticated_server_proxy # not tested
            <ServerProxy for trac.sagemath.org/login/xmlrpc>

        For convenient doctesting, this is replaced with a fake object
        during doctesting::

            sage: dev.trac._authenticated_server_proxy
            <sage.dev.trac_interface.DoctestServerProxy object at ...>
        """
        try:
            if time.time() < self.__auth_timeout:
                # default timeout is 5 minutes, like sudo
                new_timeout = time.time() + float(
                        self._config.get('password_timeout', 300))
                if new_timeout > self.__auth_timeout:
                    self.__auth_timeout = new_timeout
                return self.__authenticated_server_proxy
            else:
                del self.__authenticated_server_proxy
        except AttributeError:
            pass

        if DOCTEST_MODE:
            self.__authenticated_server_proxy = DoctestServerProxy(self)
            return self.__authenticated_server_proxy

        realm = self._config.get('realm', REALM)
        server = self._config.get('server', TRAC_SERVER_URI)

        url = urlparse.urljoin(server,
                urllib.pathname2url(os.path.join('login', 'xmlrpc')))

        while True:
            transport = DigestTransport(realm=realm, url=server,
                    username=self._username, password=self._password)
            proxy = ServerProxy(url, transport=transport)
            try:
                proxy.system.listmethods()
                break
            except urllib2.HTTPError as error:
                if error.code == 401:
                    self._UI.show("Invalid username/password pair.")
                    del self.__username
                    del self._config['password']
                else:
                    self._UI.show(
                            "Could not verify password, will try to proceed.")
                    break

        self.__authenticated_server_proxy = proxy
        return self.__authenticated_server_proxy

    @property
    def _username(self):
        """
        A lazy property to get the username on trac.

        EXAMPLES::

            sage: from sage.dev.sagedev import SageDev, Config, doctest_config
            sage: conf = doctest_config()
            sage: t = SageDev(conf).trac
            sage: t._username
            'doctest'
            sage: del conf['trac']['username']
            sage: t = SageDev(conf).trac
            sage: t._UI.append("user")
            sage: t._username
            Trac username: user
            'user'
            sage: t._username
            'user'
        """
        try:
            return self.__username
        except AttributeError:
            self.__username = self._config.get('username')
            if self.__username is None:
                self.__username = self._UI.get_input('Trac username:')
            return self.__username

    @property
    def _password(self):
        """
        A lazy property to get the username of trac.

        EXAMPLES::

            sage: from sage.dev.sagedev import SageDev, doctest_config
            sage: t = SageDev(doctest_config()).trac
            sage: t._UI.extend(["yes", "passwd", "", "pass"])
            sage: t._password
            Trac password:
            Would you like to store your password locally? (it will be stored in a plaintext file that only you may read) [yes/No/stop asking]
            'pass'
            sage: t._password
            Trac password:
            Would you like to store your password locally? (it will be stored in a plaintext file that only you may read) [yes/No/stop asking] yes
            'passwd'
            sage: t._password
            'passwd'
        """
        if self._config.get('password'):
            self.__auth_timeout = float('inf')
            return self._config['password']

        passwd = self._UI.get_password("Trac password:")

        if self._config.get('password') is None:
            r = self._UI.select("Would you like to store your password "+
                                "locally? (it will be stored in a plaintext "+
                                "file that only you may read)",
                                options=("yes","no","stop asking"), default=1)
            if r == 'yes':
                self.__auth_timeout = float('inf')
                self._config['password'] = passwd
                return passwd
            elif r == 'stop asking':
                self._config['password'] = ""

        # default timeout is 5 minutes, like sudo
        self.__auth_timeout = time.time() + float(
                self._config.get('password_timeout', 300))
        return passwd

    @property
    def sshkeys(self):
        """
        Retrieve the interface to the ssh keys stored for the user.

        EXAMPLES::

            sage: sshkeys = dev.trac.sshkeys
            sage: type(sshkeys)
            <class 'sage.dev.trac_interface.SshKeys'>
            sage: sshkeys.listusers()
            []
            sage: sshkeys.getkeys()
            []
            sage: sshkeys.setkeys(["foo", "bar"])
            0
            sage: sshkeys.getkeys()
            ['foo', 'bar']
        """
        return self._authenticated_server_proxy.sshkeys

    def create_ticket(self, summary, description, attributes={}, notify=False):
        """
        Create a ticket on trac and return the new ticket number.

        EXAMPLES::

            sage: dev.trac.create_ticket("Summary", "Description",
            ....:         {'type':'defect', 'component':'algebra'})
            14366
        """
        return self._authenticated_server_proxy.ticket.create(summary,
                description, attributes, notify)

    def edit_ticket(self, ticketnum):
        """
        edit a ticket on trac

        EXAMPLES::

            sage: import os
            sage: os.environ['EDITOR'] = 'cat'
            sage: dev.trac.edit_ticket(1000) # optional: online
            Summary: Sage does not have 10000 users yet.
            Cc:
            Type: defect
            Milestone: sage-2.10
            Component: distribution
            Priority: major
            Keywords:
            <BLANKLINE>
            ADD DESCRIPTION
            <BLANKLINE>
            <BLANKLINE>
            # Lines starting with `#` are ignored.
            # Lines starting with `Field: ` correspond to fields of
            # the trac ticket, and can be followed by text on the same line.
            # They will be assigned to the corresponding field on the trac
            # ticket.
            #
            # Lines not following this format will be put into the ticket
            # description. Trac markup is supported.
            #
            # An empty file aborts ticket creation/editing.
            Modified ticket #1000.
            sage: os.environ['EDITOR'] = 'echo "more description" >>'
            sage: dev.trac.edit_ticket(13147) # optional: online
            Modified ticket #13147.
        """
        ticketnum = int(ticketnum)
        attributes = self._get_attributes(ticketnum)

        summary = attributes.get('summary', 'No Summary')
        description = attributes.get('description', 'No Description')

        ret = self._edit_ticket_interactive(summary, description, attributes)
        if ret is None:
            return

        attributes['summary'] = ret[0]
        attributes['description'] = ret[1]
        attributes.update(ret[2])

        ticket = self._authenticated_server_proxy.ticket.update(ticketnum,
                "", attributes)[0]
        self._UI.show("Modified ticket #%s."%ticket)

    def _edit_ticket_interactive(self, summary, description, attributes):
        r"""
        edit a ticket interactively

        INPUT:

        - ``summary`` -- summary of ticket

        - ``description`` -- description of ticket

        - ``attributes`` -- dictionary containing field, value pairs

        OUTPUT: updated version of input after user has edited the ticket
        or ``None`` if editing is aborted

        TESTS::

            sage: import os
            sage: os.environ['EDITOR'] = 'cat'
            sage: tup = ('a summary', 'a description', {'Branch': 'a branch'})
            sage: tup = dev.trac._edit_ticket_interactive(*tup)
            Summary: a summary
            Branch: a branch
            a description
            <BLANKLINE>
            # Lines starting with `#` are ignored.
            # Lines starting with `Field: ` correspond to fields of
            # the trac ticket, and can be followed by text on the same line.
            # They will be assigned to the corresponding field on the trac
            # ticket.
            #
            # Lines not following this format will be put into the ticket
            # description. Trac markup is supported.
            #
            # An empty file aborts ticket creation/editing.
            sage: print tup
            ('a summary', 'a description\n', {'branch': 'a branch'})
            sage: os.environ['EDITOR'] = 'echo "added more description" >>'
            sage: print dev.trac._edit_ticket_interactive(*tup)
            ('a summary', 'a description\n\n\nadded more description\n', {'branch': 'a branch'})
            sage: os.environ['EDITOR'] = r"sed -i 's+^\(Summary: \).*$+\1new summary+'"
            sage: print dev.trac._edit_ticket_interactive(*tup)
            ('new summary', 'a description\n', {'branch': 'a branch'})
            sage: os.environ['EDITOR'] = r"sed -i 's+^\(Branch: \).*$+\1new branch+'"
            sage: print dev.trac._edit_ticket_interactive(*tup)
            ('a summary', 'a description\n', {'branch': 'new branch'})
        """
        filename = tempfile.mkstemp()[1]
        with open(filename, "w") as F:
            F.write("Summary: %s\n"%summary)
            for k,v in attributes.items():
                k = ALLOWED_FIELDS.get(k.lower())
                if k is not None:
                    F.write("%s: %s\n"%(k,v))

            if description is None or not description.strip():
                description = "\nADD DESCRIPTION\n"
            F.write(description + "\n")
            F.write(TICKET_FILE_GUIDE)

        while True:
            if self._UI.edit(filename):
                raise RuntimeError("editor exited with non-zero exit code")

            try:
                ret = _parse_ticket_file(filename)
                break
            except TicketSyntaxError as error:
                pass

            self._UI.show("TicketSyntaxError: "+error.message)

            if not self._UI.confirm("Do you want to try to fix your ticket "+
                                    "file?"):
                ret = None
                break

        os.unlink(F.name)
        return ret

    def create_ticket_interactive(self):
        """
        Interactive version of :meth:`create_ticket`.

        EXAMPLE::

            sage: from sage.dev.sagedev import SageDev, doctest_config
            sage: import os
            sage: t = SageDev(doctest_config()).trac
            sage: os.environ['EDITOR'] = 'cat'
            sage: t._UI.extend(["no"]*3)
            sage: t.create_ticket_interactive()
            Summary:
            Priority: major
            Component: PLEASE CHANGE
            Type: PLEASE CHANGE
            <BLANKLINE>
            ADD DESCRIPTION
            <BLANKLINE>
            <BLANKLINE>
            # Lines starting with `#` are ignored.
            # Lines starting with `Field: ` correspond to fields of
            # the trac ticket, and can be followed by text on the same line.
            # They will be assigned to the corresponding field on the trac
            # ticket.
            #
            # Lines not following this format will be put into the ticket
            # description. Trac markup is supported.
            #
            # An empty file aborts ticket creation/editing.
            TicketSyntaxError: no valid summary found
            Do you want to try to fix your ticket file? [Yes/no] no
            sage: os.environ['EDITOR'] = 'echo "Summary: Foo" >'
            sage: t.create_ticket_interactive()
            TicketSyntaxError: no description found
            Do you want to try to fix your ticket file? [Yes/no] no
            sage: os.environ['EDITOR'] = 'echo "Summary: Foo\nFoo: Foo\nFoo" >'
            sage: t.create_ticket_interactive()
            TicketSyntaxError: line 2: field `Foo` not supported
            Do you want to try to fix your ticket file? [Yes/no] no
            sage: os.environ['EDITOR'] = 'echo "Summary: Foo\nCc: Foo\nFoo" >'
            sage: t.create_ticket_interactive()
            Created ticket #14366.
            14366
        """
        attributes = {
                "Type":         "PLEASE CHANGE",
                "Priority":     "major",
                "Component":    "PLEASE CHANGE",
                }

        ret = self._edit_ticket_interactive("", None, attributes)

        if ret is None:
            return

        ticket = self.create_ticket(*ret)
        self._UI.show("Created ticket #%s."%ticket)

        return ticket

    def set_dependencies(self, ticket, dependencies):
        """
        Overwrites the dependencies for the given ticket.

        INPUT:

        - ``ticket`` -- an int

        - ``dependencies`` -- a list of ints
        """
        ticket = int(ticket)
        if len(dependencies) == 0:
            dep = ''
        else:
            dep = '#' + ', #'.join([str(d) for d in dependencies])
        tid, time0, time1, attributes = self._anonymous_server_proxy.ticket.get(ticket)
        olddep = attributes.get('dependencies', '')
        if dep != olddep:
            self._authenticated_server_proxy.ticket.update(tid, 'Set by SageDev: dependencies changed from %s to %s'%(olddep, dep), {'dependencies':dep})
            self._UI.show("Dependencies updated")

        # makes the trac ticket for new_ticket depend on the old_ticket
        raise NotImplementedError

    def _get_attributes(self, ticketnum):
        """
        Retrieve the properties of ticket ``ticketnum``.

        EXAMPLES::

            sage: dev.trac._get_attributes(1000) # optional: online
            {'status': 'closed',
             'changetime': <DateTime '20080110T08:28:40' at ...>,
             'description': '',
             'reporter': 'was',
             'cc': '',
             'type': 'defect',
             'milestone': 'sage-2.10',
             '_ts': '1199953720000000',
             'component': 'distribution',
             'summary': 'Sage does not have 10000 users yet.',
             'priority': 'major',
             'owner': 'was',
             'time': <DateTime '20071025T16:48:05' at ...>,
             'keywords': '',
             'resolution': 'fixed'}
        """
        return self._anonymous_server_proxy.ticket.get(int(ticketnum))[3]

    def dependencies(self, ticketnum, all=False, _seen=None):
        """
        Retrieve the dependencies of ticket ``ticketnum``.

        INPUT:

        - ``ticketnum`` -- an integer, the number of a ticket

        - ``all`` -- a boolean (default: ``False``), whether to get indirect
          dependencies of ``ticketnum``

        - ``_seen`` -- (default: ``None``), used internally in recursive calls

        EXAMPLES::

            sage: dev.trac.dependencies(1000) # optional: online (an old ticket with no dependency field)
            []
            sage: dev.trac.dependencies(13147) # optional: online
            [13579, 13681]
            sage: dev.trac.dependencies(13147,all=True) # long time, optional: online
            [13579, 13681, 13631]

        """
        # returns the list of all ticket dependencies, sorted by ticket number
        if _seen is None:
            ticketnum = int(ticketnum)
            seen = []
        elif ticketnum in _seen:
            return
        else:
            seen = _seen
        seen.append(ticketnum)
        data = self._get_attributes(ticketnum)
        if 'dependencies' not in data: return []
        dependencies = data['dependencies']
        if dependencies.strip() == '': return []
        dependencies = [a.strip(" ,;+-\nabcdefghijklmnopqrstuvwxyz") for a in data['dependencies'].split('#')]
        dependencies = [a for a in dependencies if a]
        dependencies = [int(a) if a.isdigit() else a for a in dependencies]
        if not all:
            return dependencies
        for a in dependencies:
            if isinstance(a, int):
                self.dependencies(a, True, seen)
            else:
                seen.append(a)
        if _seen is None:
            return seen[1:]

    def attachment_names(self, ticketnum):
        """
        Retrieve the names of the attachments for ticket ``ticketnum``.

        EXAMPLES::

            sage: dev.trac.attachment_names(1000) # optional: online
            ()
            sage: dev.trac.attachment_names(13147) # optional: online
            ('13147_move.patch',
             '13147_lazy.patch',
             '13147_lazy_spkg.patch',
             '13147_new.patch',
             '13147_over_13579.patch',
             'trac_13147-ref.patch',
             'trac_13147-rebased-to-13681.patch',
             'trac_13681_root.patch')
        """
        ticketnum = int(ticketnum)
        return tuple(a[0] for a in
                self._anonymous_server_proxy.ticket.listAttachments(ticketnum))

    def _set_branch(self, ticketnum, remote_branch, commit_id):
        tid = self._anonymous_server_proxy.ticket.get(int(ticketnum))[0]
        self._authenticated_server_proxy.ticket.update(tid,
                'Set by SageDev: commit %s'%(commit_id),
                {'branch':remote_branch})
