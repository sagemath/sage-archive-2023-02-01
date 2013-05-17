"""
Trac Interface
"""
import os
import re
import tempfile
import time
import urllib2

from xmlrpclib import Transport, ServerProxy

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
# An empty file aborts ticket creation.
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
        sage: _parse_ticket_file(tmp)
        ('a summary', 'some description\nmore description', {'branch': 'a branch'})
        sage: os.unlink(tmp)
    """
    lines = list(open(filename).read().splitlines())

    if all(l.rstrip().startswith('#') for l in lines if l.rstrip()):
        return

    fields = {}
    description = []

    for i, line in enumerate(lines):
        if line.startswith('#'):
            continue
        line = line.rstrip()
        i += 1 # line numbers should be indexed from 1

        m = FIELD_REGEX.match(line)
        if m and not line.startswith("sage: "):
            field = m.groups()[0]
            if not (field.lower() == 'summary' or
                    field.lower() in ALLOWED_FIELDS):
                raise TicketSyntaxError("line %s: "%i +
                                        "field `%s` not supported"%field)
            elif field.lower() in fields:
                raise TicketSyntaxError("line %s: "%i +
                                        "only one value for %s allowed"%field)
            else:
                fields[field.lower()] = m.groups()[1].strip()
                continue

        if line != "[Description]":
            description.append(line)

    # no syntax errors in file
    try:
        summary = fields.pop('summary')
    except KeyError:
        summary = None

    if not summary:
        raise TicketSyntaxError("no valid summary found")
    elif not "".join(description):
        raise TicketSyntaxError("no description found")
    else:
        return summary, "\n".join(description), fields

class DigestTransport(object, Transport):
    """
    Handles an HTTP transaction to an XML-RPC server.

    EXAMPLES::

        sage: from sage.env import REALM, TRAC_SERVER_URI
        sage: sage.dev.trac_interface.DigestTransport(REALM, TRAC_SERVER_URI+"/xmlrpc")
        <sage.dev.trac_interface.DigestTransport object at ...>
    """
    def __init__(self, realm, url, username=None, password=None, **kwds):
        """
        Initialization.

        EXAMPLES::

            sage: from sage.env import REALM, TRAC_SERVER_URI
            sage: type(sage.dev.trac_interface.DigestTransport(REALM, TRAC_SERVER_URI+"/xmlrpc"))
            <class 'sage.dev.trac_interface.DigestTransport'>
        """
        Transport.__init__(self, **kwds)

        authhandler = urllib2.HTTPDigestAuthHandler()
        if username and password:
            authhandler.add_password(realm, url, username, password)

        self.opener = urllib2.build_opener(authhandler)

    def request(self, host, handler, request_body, verbose=0):
        """
        Issue an XML-RPC request.

        EXAMPLES::

            sage: from sage.env import REALM, TRAC_SERVER_URI
            sage: d = sage.dev.trac_interface.DigestTransport(REALM, TRAC_SERVER_URI+"/xmlrpc")
            sage: d.request # not tested
        """
        self.verbose = verbose

        headers = {'Content-type': 'text/xml'}
        data = request_body
        req = urllib2.Request('http://' + host + handler, data, headers)

        response = self.opener.open(req)

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
        if 'trac' not in sagedev._config:
            sagedev._config['trac'] = {}
        self._config = sagedev._config['trac']

        # Caches for the analogous single-underscore properties
        self.__anonymous_server_proxy = None
        self.__authenticated_server_proxy = None

        self.__passwd = None
        self.__passwd_timeout = None

    @property
    def _anonymous_server_proxy(self):
        """
        Lazy wrapper around a non-authenticated XML-RPC interface to trac.

        .. NOTE::

            Unlike the authenticated server proxy, this is not replaced with a
            fake proxy for doctesting. All doctests using it should therefore
            be labeled as optional ``online``

        EXAMPLES::

            sage: dev.trac._anonymous_server_proxy
            <ServerProxy for trac.tangentspace.org/sage_trac/xmlrpc>
        """
        if self.__anonymous_server_proxy is None:
            realm = REALM
            if "realm" in self._config:
                realm = self._config["realm"]
            server = TRAC_SERVER_URI
            if "server" in self._config:
                server = self._config["server"]
            if server[-1] != '/': server += '/'

            transport = DigestTransport(realm, server)
            self.__anonymous_server_proxy = ServerProxy(server + 'xmlrpc', transport=transport)
        return self.__anonymous_server_proxy

    @property
    def _authenticated_server_proxy(self):
        """
        Get an XML-RPC proxy object that is authenticated using the users
        username and password.

        EXAMPLES::

            sage: dev.trac._authenticated_server_proxy # not tested
            <ServerProxy for trac.tangentspace.org/sage_trac/login/xmlrpc>

        For convenient doctesting, this is replaced with a fake object for the user ``'doctest'``::

            sage: dev.trac._authenticated_server_proxy
            <sage.dev.trac_interface.DoctestServerProxy object at ...>
        """
        config = self._config

        if self.__authenticated_server_proxy is None:
            realm = REALM
            if "realm" in self._config:
                realm = self._config["realm"]
            server = TRAC_SERVER_URI
            if "server" in self._config:
                server = self._config["server"]
            if server[-1] != '/': server += '/'

            username = self._username
            if username == "doctest":
                return DoctestServerProxy(self)
            else:
                transport = DigestTransport(realm, server, username, self._password)
                self.__authenticated_server_proxy = ServerProxy(server + 'login/xmlrpc', transport=transport)

        return self.__authenticated_server_proxy


    @property
    def _username(self):
        """
        A lazy property to get the username on trac.

        EXAMPLES::

            sage: from sage.dev.sagedev import SageDev, Config, doctest_config
            sage: conf = doctest_config()
            sage: del conf['trac']
            sage: t = SageDev(conf).trac
            sage: t._UI.append("user")
            sage: t._username
            Please enter your trac username: user
            'user'
            sage: t._username
            'user'
            sage: SageDev(Config(conf._devrc)).trac._username
            'user'
        """
        if 'username' not in self._config:
            self._config['username'] = self._UI.get_input("Please enter your trac username:")
        return self._config['username']

    @property
    def _password(self):
        """
        A lazy property to get the username of trac.

        EXAMPLES::

            sage: from sage.dev.sagedev import SageDev, doctest_config
            sage: t = SageDev(doctest_config()).trac
            sage: t._UI.extend(["yes", "passwd", "passwd", "", "pass", "pass"])
            sage: t._password
            Please enter your trac password:
            Please confirm your trac password:
            Do you want your password to be stored on your local system? (your password will be stored in plaintext in a file only readable by you) [yes/No/stop asking]
            'pass'
            sage: t._password
            'pass'
            sage: import time      # long time
            sage: time.sleep(1)    # long time
            sage: t._password      # long time
            Please enter your trac password:
            Please confirm your trac password:
            Do you want your password to be stored on your local system? (your password will be stored in plaintext in a file only readable by you) [yes/No/stop asking] yes
            'passwd'
            sage: time.sleep(1)    # long time
            sage: t._password      # long time
            'passwd'
        """
        if self._config.get('password'):
            return self._config['password']

        if self.__passwd_timeout is not None:
            if time.time() < self.__passwd_timeout:
                return self.__passwd

        while True:
            passwd = self._UI.get_password("Please enter your trac password:")
            if (self._UI.get_password("Please confirm your trac password:")
                    == passwd):
                break
            else:
                self._UI.show("Passwords do not agree.")

        self.__passwd = passwd
        self.__passwd_timeout = time.time()
        # default timeout is 15 minutes, like sudo
        self.__passwd_timeout += float(
                self._config.get('password_timeout', 900))

        if self._config.get('password') is None:
            r = self._UI.select("Do you want your password to be stored on "+
                                "your local system? (your password will be "+
                                "stored in plaintext in a file only readable "+
                                "by you)",
                                options=("yes","no","stop asking"), default=1)
            if r == 'yes':
                self._config['password'] = passwd
            elif r == 'stop asking':
                self._config['password'] = ""

        return self.__passwd

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
            sage: sshkeys.setkeys(["foo","bar"])
            0
            sage: sshkeys.getkeys()
            ['foo', 'bar']
        """
        return self._authenticated_server_proxy.sshkeys

    def create_ticket(self, summary, description, attributes={}, notify=False):
        """
        Create a ticket on trac and return the new ticket number.

        EXAMPLES::

            sage: dev.trac.create_ticket("Summary","Description",{'type':'defect','component':'algebra'})
            14366
        """
        return self._authenticated_server_proxy.ticket.create(summary, description, attributes, notify)

    def edit_ticket(self, ticketnum):
        attributes = self._get_attributes(ticketnum)

        summary = "No Summary"
        if 'summary' in attributes:
            summary = attributes['summary']
        summary += "(can not be changed)"

        description = "No Description"
        if 'description' in attributes:
            description = attributes['description']

        while True:
            try:
                x = self._edit_ticket_interactive(summary, description, attributes)
                if x is None: return
                summary, description, attributes = x
                attributes['description'] = description
                self._authenticated_server_proxy.ticket.update(ticketnum, "", attributes)
            except StandardError:
                self._UI.show("Ticket editing failed: %s"%e)
                if self._UI.confirm("Do you want to try to fix your ticket file?"): continue
                else: return None

    def _edit_ticket_interactive(self, summary, description, attributes):
        """
        edit a ticket interactively
        """
        filename = tempfile.mkstemp()[1]
        with open(filename, "w") as F:
            F.write("Summary: %s\n"%summary)
            for k,v in attributes.items():
                k = ALLOWED_FIELDS.get(k.lower())
                if k is not None:
                    F.write("%s: %s\n"%(k,v))

            if description is None or not description.strip():
                description = "\n[Description]\n"
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
            sage: t._UI.extend(["yes"]+["no", "yes"]*3)
            sage: t.create_ticket_interactive()
            Do you want to create a new ticket? [Yes/no] yes
            Summary:
            Priority: major
            Keywords:
            Type: defect
            <BLANKLINE>
            [Description]
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
            # An empty file aborts ticket creation.
            TicketSyntaxError: no valid summary found
            Do you want to try to fix your ticket file? [Yes/no] no
            sage: os.environ['EDITOR'] = 'echo "Summary: Foo" >'
            sage: t.create_ticket_interactive()
            Do you want to create a new ticket? [Yes/no] yes
            TicketSyntaxError: no description found
            Do you want to try to fix your ticket file? [Yes/no] no
            sage: os.environ['EDITOR'] = 'echo "Summary: Foo\nFoo\nFoo: Foo" >'
            sage: t.create_ticket_interactive()
            Do you want to create a new ticket? [Yes/no] yes
            TicketSyntaxError: line 3: field `Foo` not supported
            Do you want to try to fix your ticket file? [Yes/no] no
            sage: os.environ['EDITOR'] = 'echo "Summary: Foo\nFoo\nCc: Foo" >'
            sage: t.create_ticket_interactive()
            Do you want to create a new ticket? [Yes/no] yes
            Created ticket #14366.
            14366
        """
        if self._UI.confirm("Do you want to create a new ticket?"):
            summary, description, attributes = "","\n",{"Type":"defect","Priority":"major","Keywords":""}
            while True:
                try:
                    ret = self._edit_ticket_interactive(summary, description, attributes)
                    if not ret: return None
                    summary, description, attributes = ret
                    ticket = self.create_ticket(summary, description, attributes)
                    self._UI.show("Created ticket #%s."%ticket)
                    return ticket
                except StandardError as e:
                        self._UI.show("Ticket creation failed: %s"%e)
                        if self._UI.confirm("Do you want to try to fix your ticket file?"): continue
                        else: return None
        assert(False)

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
            {'_ts': '1199953720000000',
             'cc': '',
             'changetime': <DateTime '20080110T08:28:40' at ...>,
             'component': 'distribution',
             'description': '',
             'keywords': '',
             'milestone': 'sage-2.10',
             'owner': 'was',
             'priority': 'major',
             'reporter': 'was',
             'resolution': 'fixed',
             'status': 'closed',
             'summary': 'Sage does not have 10000 users yet.',
             'time': <DateTime '20071025T16:48:05' at ...>,
             'type': 'defect'}

        """
        ticketnum = int(ticketnum)
        return self._anonymous_server_proxy.ticket.get(ticketnum)[3]

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
            []
            sage: dev.trac.attachment_names(13147) # optional: online
            ['13147_move.patch', '13147_lazy.patch', '13147_lazy_spkg.patch', '13147_new.patch', '13147_over_13579.patch', 'trac_13147-ref.patch', 'trac_13147-rebased-to-13681.patch', 'trac_13681_root.patch']
        """
        ticketnum = int(ticketnum)
        attachments = self._anonymous_server_proxy.ticket.listAttachments(ticketnum)
        return [a[0] for a in attachments]

    def _set_branch(self, ticketnum, remote_branch, commit_id):
        ticketnum = int(ticketnum)
        tid, time0, time1, attributes = self._anonymous_server_proxy.ticket.get(ticketnum)
        self._authenticated_server_proxy.ticket.update(tid, 'Set by SageDev: commit %s'%(commit_id), {'branch':remote_branch})
