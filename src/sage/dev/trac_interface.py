r"""
Trac Interface

This module provides an interface to access sage's issue tracker 'trac' through
its RPC interface.

AUTHORS:

- David Roe, Julian Rueth, R. Andrew Ohana, Robert Bradshaw, Timo Kluck:
  initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 David Roe <roed.math@gmail.com>
#                          Julian Rueth <julian.rueth@fsfe.org>
#                          R. Andrew Ohana <andrew.ohana@gmail.com>
#                          Robert Bradshaw <robertwb@gmail.com>
#                          Timo Kluck <tkluck@infty.nl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import re
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
        "status":           "Status",
        "type":             "Type",
        "work issues":      "Work issues",
        }
TICKET_FILE_GUIDE = r"""
# Lines starting with `#` are ignored.
# Lines at the beginning of this file starting with `Field: ` correspond to
# fields of the trac ticket, and can be followed by text on the same line.
# They will be assigned to the corresponding field on the trac ticket.
#
# Lines not following this format will be put into the ticket description. Trac
# markup is supported.
#
# An empty file aborts ticket creation/editing.
"""
COMMENT_FILE_GUIDE = r"""
# Lines starting with `#` are ignored.
# An empty file aborts the comment.
"""

class TicketSyntaxError(SyntaxError): # we don't want to catch normal syntax errors
    r"""
    A syntax error when parsing a ticket description modified by the user.

    EXAMPLES::

        sage: from sage.dev.trac_interface import TicketSyntaxError
        sage: raise TicketSyntaxError()
        Traceback (most recent call last):
        ...
        TicketSyntaxError: None

    """

class TracInterface(object):
    r"""
    Wrapper around the XML-RPC interface of trac.

    EXAMPLES::

        sage: from sage.dev.test.config import DoctestConfig
        sage: from sage.dev.test.user_interface import DoctestUserInterface
        sage: from sage.dev.trac_interface import TracInterface
        sage: config = DoctestConfig()
        sage: trac = TracInterface(config['trac'], DoctestUserInterface(config['UI']))
        sage: trac
        <sage.dev.trac_interface.TracInterface object at 0x...>

    """
    def __init__(self, config, UI):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.trac_interface import TracInterface
            sage: config = DoctestConfig()
            sage: trac = TracInterface(config['trac'], DoctestUserInterface(config['UI']))
            sage: type(trac)
            <class 'sage.dev.trac_interface.TracInterface'>

        """
        self._UI = UI
        self._config = config

        self.__anonymous_server_proxy = None
        self.__authenticated_server_proxy = None

        self.__username = None
        self.__password = None
        self.__auth_timeout = None

    @property
    def _username(self):
        r"""
        A lazy property to get the username on trac.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.trac_interface import TracInterface
            sage: config = DoctestConfig()
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = TracInterface(config['trac'], UI)
            sage: trac._username # username is read from config
            'doctest'
            sage: trac.reset_username()
            sage: UI.append('doctest2')
            sage: trac._username # user is prompted for a username
            Trac username: doctest2
            'doctest2'
            sage: config['trac']['username']
            'doctest2'

        """
        if self.__username is None:
            self.__username = self._config.get('username', None)

        if self.__username is None:
            self.__username = self._config['username'] = self._UI.get_input('Trac username:')
            from user_interface import INFO
            self._UI.show("Your trac username has been written to a configuration file for future sessions. To reset your username, use `dev.trac.reset_username()`.", INFO)

        return self.__username

    def reset_username(self):
        r"""
        Reset username and password stored in this object and in the
        configuration.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.trac_interface import TracInterface
            sage: config = DoctestConfig()
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = TracInterface(config['trac'], UI)
            sage: trac.reset_username()
            sage: UI.append("doctest2")
            sage: trac._username
            Trac username: doctest2
            'doctest2'

        """
        self.__username = None
        if 'username' in self._config:
            del self._config['username']

        self.reset_password()

    @property
    def _password(self):
        r"""
        A lazy property to get the password for trac.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.trac_interface import TracInterface
            sage: config = DoctestConfig()
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = TracInterface(config['trac'], UI)
            sage: UI.append('')
            sage: UI.append('secret')
            sage: trac._password
            Trac password:
            Should I store your password in a configuration file for future sessions? (This configuration file might be readable by privileged users on this system.) [yes/No]
            'secret'
            sage: trac._password # password is stored for some time, so there is no need to type it immediately afterwards
            'secret'
            sage: config['trac']['password']
            Traceback (most recent call last):
            ...
            KeyError: 'password'
            sage: trac.reset_password()

            sage: UI.append('y')
            sage: UI.append('secret')
            sage: trac._password
            Trac password:
            Should I store your password in a configuration file for future sessions? (This configuration file might be readable by privileged users on this system.) [yes/No] y
            'secret'
            sage: config['trac']['password']
            'secret'
            sage: trac._password
            'secret'

        """
        self._check_password_timeout()

        if self.__password is None:
            self.__password = self._config.get('password', None)

        if self.__password is None:
            self.__password = self._UI.get_password('Trac password:')
            store_password = self._config.get('store_password', None)
            if store_password is None:
                store_password = "yes" if self._UI.confirm("Should I store your password in a configuration file for future sessions? (This configuration file might be readable by privileged users on this system.)", default=False) else "no"
                if store_password == "no":
                    self._config['store_password'] = store_password # remember the user's decision (if negative) and do not ask every time

            if store_password == "yes":
                self._config['password'] = self.__password
                from sage.dev.user_interface import INFO
                self._UI.show("Your trac password has been written to a configuration file. To reset your password, use `dev.trac.reset_password()`.", INFO)

        self._postpone_password_timeout()
        return self.__password

    def reset_password(self):
        r"""
        Reset password stored in this object and in the configuration.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.trac_interface import TracInterface
            sage: config = DoctestConfig()
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = TracInterface(config['trac'], UI)
            sage: UI.append('y')
            sage: UI.append('secret')
            sage: trac._password
            Trac password:
            Should I store your password in a configuration file for future sessions? (This configuration file might be readable by privileged users on this system.) [yes/No] y
            'secret'
            sage: config['trac']['password']
            'secret'
            sage: trac.reset_password()
            sage: config['trac']['password']
            Traceback (most recent call last):
            ...
            KeyError: 'password'

        """
        self.__password = None
        self.__authenticated_server_proxy = None
        self.__auth_timeout = None
        if 'password' in self._config:
            del self._config['password']
        if 'store_password' in self._config:
            del self._config['store_password']

    def _check_password_timeout(self):
        r"""
        Reset all attributes that depend on the saved password if it has timed
        out (usually after 5 minutes without using it).

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.trac_interface import TracInterface
            sage: config = DoctestConfig()
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = TracInterface(config['trac'], UI)
            sage: UI.append('')
            sage: UI.append('secret')
            sage: config['trac']['password_timeout'] = 0
            sage: trac._password
            Trac password:
            Should I store your password in a configuration file for future sessions? (This configuration file might be readable by privileged users on this system.) [yes/No]
            'secret'
            sage: UI.append('secret')
            sage: trac._password # indirect doctest
            Trac password:
            'secret'
            sage: trac.reset_password()
            sage: UI.append('y')
            sage: UI.append('secret')
            sage: trac._password
            Trac password:
            Should I store your password in a configuration file for future sessions? (This configuration file might be readable by privileged users on this system.) [yes/No] y
            'secret'
            sage: trac._password # the timeout has no effect if the password can be read from the configuration file
            'secret'

        """
        import time
        if self.__auth_timeout is None or time.time() >= self.__auth_timeout:
            self.__password = None
            self.__authenticated_server_proxy = None
            self.__auth_timeout = None

    def _postpone_password_timeout(self):
        r"""
        Postpone the password timeout.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.trac_interface import TracInterface
            sage: config = DoctestConfig()
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = TracInterface(config['trac'], UI)
            sage: UI.append('')
            sage: UI.append('secret')
            sage: trac._password
            Trac password:
            Should I store your password in a configuration file for future sessions? (This configuration file might be readable by privileged users on this system.) [yes/No]
            'secret'
            sage: trac._password # indirect doctest
            'secret'

        """
        import time
        new_timeout = time.time() + float(self._config.get('password_timeout', 300))

        if self.__auth_timeout is None or new_timeout > self.__auth_timeout:
            self.__auth_timeout = new_timeout

    @property
    def _anonymous_server_proxy(self):
        """
        Return a non-authenticated XML-RPC interface to trac.

        .. NOTE::

            Unlike the authenticated server proxy, this can be used in
            doctesting. However, all doctests relying on it talking to the
            actual trac server should be marked as ``optional: internet``.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.trac_interface import TracInterface
            sage: config = DoctestConfig()
            sage: trac = TracInterface(config['trac'], DoctestUserInterface(config['UI']))
            sage: trac._anonymous_server_proxy
            <ServerProxy for trac.sagemath.org/xmlrpc>

        """
        if self.__anonymous_server_proxy is None:
            from sage.env import TRAC_SERVER_URI
            server = self._config.get('server', TRAC_SERVER_URI)
            import urlparse
            url = urlparse.urljoin(server, 'xmlrpc')
            from digest_transport import DigestTransport
            transport = DigestTransport()
            from xmlrpclib import ServerProxy
            self.__anonymous_server_proxy = ServerProxy(url, transport=transport)

        return self.__anonymous_server_proxy

    @property
    def _authenticated_server_proxy(self):
        r"""
        Get an XML-RPC proxy object that is authenticated using the users
        username and password.

        .. NOTE::

            To make sure that doctests do not tamper with the live trac server,
            it is an error to access this property during a doctest.

        EXAMPLES::

            sage: dev.trac._authenticated_server_proxy # not tested
            Trac username: username
            Trac password:
            Should I store your password in a configuration file for future sessions? (This configuration file might be readable by privileged users on this system.) [yes/No]
            <ServerProxy for trac.sagemath.org/login/xmlrpc>

        TESTS::

            sage: dev.trac._authenticated_server_proxy
            Traceback (most recent call last):
            ...
            AssertionError: doctest tried to access an authenticated session to trac

        """
        import sage.doctest
        assert not sage.doctest.DOCTEST_MODE, "doctest tried to access an authenticated session to trac"

        self._check_password_timeout()

        if self.__authenticated_server_proxy is None:
            from sage.env import REALM
            realm = self._config.get('realm', REALM)
            from sage.env import TRAC_SERVER_URI
            server = self._config.get('server', TRAC_SERVER_URI)

            import os, urllib, urllib2, urlparse
            url = urlparse.urljoin(server, urllib.pathname2url(os.path.join('login', 'xmlrpc')))

            while True:
                from xmlrpclib import ServerProxy
                from digest_transport import DigestTransport
                transport = DigestTransport()
                transport.add_authentication(realm=realm, url=server, username=self._username, password=self._password)
                proxy = ServerProxy(url, transport=transport)
                try:
                    proxy.system.listMethods()
                    break
                except urllib2.HTTPError as error:
                    if error.code == 401:
                        self._UI.show("Invalid username/password")
                        self.reset_username()
                    else:
                        self._UI.show("Could not verify password, will try to proceed.")
                        break

            self.__authenticated_server_proxy = proxy

        self._postpone_password_timeout()

        return self.__authenticated_server_proxy

    def create_ticket(self, summary, description, attributes={}):
        r"""
        Create a ticket on trac and return the new ticket number.

        .. SEEALSO::

            :meth:`create_ticket_interactive`

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: trac.create_ticket('Summary', 'Description', {'type':'defect', 'component':'algebra'})
            1

        """
        return self._authenticated_server_proxy.ticket.create(summary, description, attributes)

    def add_comment(self, ticket, comment):
        r"""
        Add ``comment`` to ``ticket`` on trac.

        .. SEEALSO::

            :meth:`add_comment_interactive`

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: ticket = trac.create_ticket('Summary', 'Description', {'type':'defect', 'component':'algebra'})
            sage: trac.add_comment(ticket, "a comment")

        """
        ticket = int(ticket)
        attributes = self._get_attributes(ticket)
        self._authenticated_server_proxy.ticket.update(ticket, comment, attributes)

    def _get_attributes(self, ticket):
        r"""
        Retrieve the properties of ``ticket``.

        EXAMPLES::

            sage: dev.trac._get_attributes(1000) # optional: internet
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
        return self._anonymous_server_proxy.ticket.get(int(ticket))[3]

    def _branch_for_ticket(self, ticket):
        r"""
        Return the branch field for ``ticket`` or ``None`` if it is not set.

        INPUT:

        - ``ticket`` -- an int

        EXAMPLES::

            sage: dev.trac._branch_for_ticket(1000) is None # optional: internet
            True

        """
        attributes = self._get_attributes(ticket)
        if 'branch' in attributes:
            return attributes['branch'] or None
        else:
            return None

    def dependencies(self, ticket, recurse=False, seen=None):
        r"""
        Retrieve dependencies of ``ticket``, sorted by ticket number.

        INPUT:

        - ``ticket`` -- an integer, the number of the ticket

        - ``recurse`` -- a boolean (default: ``False``), whether to get
          indirect dependencies of ``ticket``

        - ``seen`` -- a list (default: ``[]``), used internally to implement
          ``recurse``

       EXAMPLES::

            sage: dev.trac.dependencies(1000)            # optional: internet (an old ticket with no dependency field)
            []
            sage: dev.trac.dependencies(13147)           # optional: internet
            [13579, 13681]
            sage: dev.trac.dependencies(13147, recurse=True) # long time, optional: internet
            [13579, 13631, 13681]

        """
        ticket = int(ticket)

        if seen is None:
            seen = []

        if ticket in seen:
            return []

        seen.append(ticket)
        dependencies = self._get_attributes(ticket).get('dependencies','').strip()
        dependencies = dependencies.split(',')
        dependencies = [dep.strip() for dep in dependencies]
        dependencies = [dep for dep in dependencies if dep]
        if not all(dep[0]=="#" for dep in dependencies):
            raise RuntimeError("malformatted dependency on ticket `%s`"%ticket)
        dependencies = [dep[1:] for dep in dependencies]
        try:
            dependencies = [int(dep) for dep in dependencies]
        except ValueError:
            raise RuntimeError("malformatted dependency on ticket `%s`"%ticket)

        if recurse:
            for dep in dependencies:
                self.dependencies(dep, recurse, seen)
        else:
            seen.extend(dependencies)

        ret = sorted(seen)
        ret.remove(ticket)
        return ret

    def attachment_names(self, ticket):
        """
        Retrieve the names of the attachments for ``ticket``.

        EXAMPLES::

            sage: dev.trac.attachment_names(1000) # optional: internet
            ()
            sage: dev.trac.attachment_names(13147) # optional: internet
            ('13147_move.patch',
             '13147_lazy.patch',
             '13147_lazy_spkg.patch',
             '13147_new.patch',
             '13147_over_13579.patch',
             'trac_13147-ref.patch',
             'trac_13147-rebased-to-13681.patch',
             'trac_13681_root.patch')
        """
        ticket = int(ticket)
        return tuple(a[0] for a in self._anonymous_server_proxy.ticket.listAttachments(ticket))

    def add_comment_interactive(self, ticket, comment=''):
        r"""
        Add a comment to ``ticket`` on trac.

        INPUT:

        - ``comment`` -- a string (default: ``''``), the default value for the
          comment to add.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: ticket = trac.create_ticket('Summary', 'Description', {'type':'defect', 'component':'algebra'})

            sage: UI.append("# empty comment")
            sage: trac.add_comment_interactive(ticket)
            Traceback (most recent call last):
            ...
            OperationCancelledError: comment creation aborted

            sage: UI.append("a comment")
            sage: trac.add_comment_interactive(ticket)

        """
        ticket = int(ticket)

        attributes = self._get_attributes(ticket)

        import tempfile
        fd, filename = tempfile.mkstemp()
        import os
        with os.fdopen(fd, "w") as F:
            F.write(comment)
            F.write("\n")
            F.write(COMMENT_FILE_GUIDE)

        self._UI.edit(filename)

        comment = list(open(filename).read().splitlines())
        comment = [line for line in comment if not line.startswith("#")]
        if all([line.strip()=="" for line in comment]):
            from user_interface_error import OperationCancelledError
            raise OperationCancelledError("comment creation aborted")
        comment = "\n".join(comment)

        url = self._authenticated_server_proxy.ticket.update(ticket, comment, attributes)

        from user_interface import INFO
        self._UI.show("Your comment has been recorded: %s"%url, INFO)

    def edit_ticket_interactive(self, ticket):
        """
        Edit ``ticket`` on trac.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: ticket = trac.create_ticket('Summary', 'Description', {'type':'defect', 'component':'algebra'})

            sage: UI.append("# empty")
            sage: trac.edit_ticket_interactive(ticket)
            Traceback (most recent call last):
            ...
            OperationCancelledError: ticket edit aborted

            sage: UI.append("Summary: summary\ndescription\n")
            sage: trac.edit_ticket_interactive(ticket)

        """
        ticket = int(ticket)
        attributes = self._get_attributes(ticket)

        summary = attributes.get('summary', 'No Summary')
        description = attributes.get('description', 'No Description')

        ret = self._edit_ticket_interactive(summary, description, attributes)
        if ret is None:
            from user_interface_error import OperationCancelledError
            raise OperationCancelledError("edit aborted")

        attributes['summary'] = ret[0]
        attributes['description'] = ret[1]
        attributes.update(ret[2])

        url = self._authenticated_server_proxy.ticket.update(ticket, "", attributes)
        from user_interface import INFO
        self._UI.show("Ticket modified: %s"%url, INFO)

    def _edit_ticket_interactive(self, summary, description, attributes):
        r"""
        Helper method for :meth:`edit_ticket_interactive` and
        :meth:`create_ticket_interactive`.

        INPUT:

        - ``summary`` -- a string, summary of ticket

        - ``description`` -- a string, description of ticket

        - ``attributes`` -- dictionary containing field, value pairs

        OUTPUT:

        A tuple ``(summary, description, attributes)``, the updated version of
        input after user has edited the ticket.

        TESTS::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.trac_interface import TracInterface
            sage: config = DoctestConfig()
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = TracInterface(config['trac'], UI)
            sage: UI.append("# abort")
            sage: trac._edit_ticket_interactive('summary', 'description', {'branch':'branch1'})
            Traceback (most recent call last):
            ...
            OperationCancelledError: ticket edit aborted

            sage: UI.append("Summary: new summary\nBranch: branch2\nnew description")
            sage: trac._edit_ticket_interactive('summary', 'description', {'branch':'branch1'})
            ('new summary', 'new description', {'branch': 'branch2'})

            sage: UI.append("Summary: new summary\nBranch: branch2\nnew description")
            sage: UI.append("")
            sage: UI.append("Summary: new summary\nInvalid: branch2\nnew description")
            sage: trac._edit_ticket_interactive('summary', 'description', {'branch':'branch1'})
            TicketSyntaxError: line 2: field `Invalid` not supported
            Do you want to try to fix your ticket file? [Yes/no]
            ('new summary', 'new description', {'branch': 'branch2'})

        """
        import tempfile, os
        fd, filename = tempfile.mkstemp()
        try:
            with os.fdopen(fd, "w") as F:
                F.write("Summary: %s\n"%summary)
                for k,v in attributes.items():
                    k = ALLOWED_FIELDS.get(k.lower())
                    if k is not None:
                        F.write("%s: %s\n"%(k,v))

                if description is None or not description.strip():
                    description = "\nADD DESCRIPTION\n"
                F.write("\n" + description + "\n")
                F.write(TICKET_FILE_GUIDE)

            while True:
                try:
                    self._UI.edit(filename)
                    ret = self._parse_ticket_file(filename)
                    break
                except (RuntimeError, TicketSyntaxError) as error:
                    pass

                self._UI.show("TicketSyntaxError: "+error.message)

                if not self._UI.confirm("Do you want to try to fix your ticket file?", default=True):
                    ret = None
                    break

            if ret is None:
                from user_interface_error import OperationCancelledError
                raise OperationCancelledError("ticket edit aborted")

        finally:
            os.unlink(filename)

        return ret

    def create_ticket_interactive(self):
        r"""
        Drop user into an editor for creating a ticket.

        EXAMPLE::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: UI.append("Summary: summary\nType: defect\nPriority: minor\nComponent: algebra\ndescription")
            sage: trac.create_ticket_interactive()
            1

        """
        attributes = {
                "Type":         "PLEASE CHANGE",
                "Priority":     "major",
                "Component":    "PLEASE CHANGE",
                }

        ret = self._edit_ticket_interactive("", None, attributes)

        if ret is None:
            from user_interface_error import OperationCancelledError
            raise OperationCancelledError("ticket creation aborted")

        ticket = self.create_ticket(*ret)
        import urlparse
        from sage.env import TRAC_SERVER_URI
        ticket_url = urlparse.urljoin(self._config.get('server', TRAC_SERVER_URI), str(ticket))
        from user_interface import INFO
        self._UI.show("Created ticket #%s (%s)."%(ticket, ticket_url), INFO)
        return ticket

    @classmethod
    def _parse_ticket_file(cls, filename):
        r"""
        Parse ticket file ``filename``, helper for
        :meth:`create_ticket_interactive` and :meth:`edit_ticket_interactive`.

        OUTPUT:

        ``None`` if the filename contains only comments; otherwise a triple
        ``(summary, description, attributes)``, where ``summary`` is a string
        consisting of the ticket's summary, ``description`` is a string
        containing the ticket's description, and ``attributes`` is a dictionary
        with additional fields of the ticket.

        TESTS::

            sage: from sage.dev.trac_interface import TracInterface
            sage: import os, tempfile
            sage: tmp = tempfile.mkstemp()[1]
            sage: with open(tmp, 'w') as f:
            ....:     f.write("no summary\n")
            sage: TracInterface._parse_ticket_file(tmp)
            Traceback (most recent call last):
            ...
            TicketSyntaxError: no valid summary found
            sage: with open(tmp, 'w') as f:
            ....:     f.write("summary:no description\n")
            sage: TracInterface._parse_ticket_file(tmp)
            Traceback (most recent call last):
            ...
            TicketSyntaxError: no description found
            sage: with open(tmp, 'w') as f:
            ....:     f.write("summary:double summary\n")
            ....:     f.write("summary:double summary\n")
            sage: TracInterface._parse_ticket_file(tmp)
            Traceback (most recent call last):
            ...
            TicketSyntaxError: line 2: only one value for summary allowed
            sage: with open(tmp, 'w') as f:
            ....:     f.write("bad field:bad field entry\n")
            sage: TracInterface._parse_ticket_file(tmp)
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
            sage: TracInterface._parse_ticket_file(tmp)
            ('a summary', 'some description\nmore description', {'branch': 'a branch'})
            sage: with open(tmp, 'w') as f:
            ....:     f.write("summary:a summary\n")
            ....:     f.write("some description\n")
            ....:     f.write("branch:a branch\n")
            sage: TracInterface._parse_ticket_file(tmp)
            ('a summary', 'some description\nbranch:a branch', {})
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
            return summary, "\n".join(description), fields
