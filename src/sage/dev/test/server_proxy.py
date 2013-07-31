r"""
Trac Server Proxy for Doctesting

This module provides substitutes for the server proxy used by
:class:`sage.dev.trac_interface.TracInterface` for doctesting.

AUTHORS:

- TODO: add authors from github's history and trac's history

"""
#*****************************************************************************
#       Copyright (C) 2013 TODO
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.dev.trac_interface

class DoctestServerProxy(object):
    r"""
    A server proxy which can be used by
    :meth:`sage.dev.test.trac_interface.DoctestTracInterface._anonymous_server_proxy`
    for doctesting.

    EXAMPLES::

        sage: from sage.dev.test.server_proxy import DoctestServerProxy
        sage: from sage.dev.test.trac_server import DoctestTracServer
        sage: DoctestServerProxy(DoctestTracServer())
         <sage.dev.test.server_proxy.DoctestServerProxy at 0x...>

    """
    def __init__(self, server):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: type(trac._anonymous_server_proxy)
            <class 'sage.dev.test.server_proxy.DoctestServerProxy'>

        """
        self._server = server

        self.ticket = DoctestTicketProxy(self)

    def _check_authentication(self, privilege):
        r"""
        Check whether the user has sufficient permissions to perform an action
        which requires ``privilege``.

        EXAMPLES::

            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: trac._anonymous_server_proxy._check_authentication("TICKET_CREATE")
            Traceback (most recent call last):
            ...
            Fault: <Fault 403: "TICKET_CREATE privileges are required to perform this operation. You don't have the required permissions.">
            sage: trac._authenticated_server_proxy._check_authentication("TICKET_CREATE")

        """
        import xmlrpclib
        raise xmlrpclib.Fault(403, "%s privileges are required to perform this operation. You don't have the required permissions."%privilege)

class AuthenticatedDoctestServerProxy(DoctestServerProxy):
    r"""
    A server proxy which can be used by
    :meth:`sage.dev.test.trac_interface.DoctestTracInterface._anonymous_server_proxy`
    for doctesting.

    EXAMPLES::

        sage: from sage.dev.test.server_proxy import AuthenticatedDoctestServerProxy
        sage: from sage.dev.test.trac_server import DoctestTracServer
        sage: AuthenticatedDoctestServerProxy(DoctestTracServer(), 'username', 'password')
        <sage.dev.test.server_proxy.AuthenticatedDoctestServerProxy at 0x...>

    """
    def __init__(self, server, username, password):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: type(trac._authenticated_server_proxy)
            <class 'sage.dev.test.server_proxy.AuthenticatedDoctestServerProxy'>

        """
        DoctestServerProxy.__init__(self, server)

        self._username = username
        self._password = password

    def _check_authentication(self, privilege):
        r"""
        Check whether the user has sufficient permissions to perform an action
        which requires ``privilege``.

        EXAMPLES::

            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: trac._anonymous_server_proxy._check_authentication("TICKET_CREATE")
            Traceback (most recent call last):
            ...
            Fault: <Fault 403: "TICKET_CREATE privileges are required to perform this operation. You don't have the required permissions.">
            sage: trac._authenticated_server_proxy._check_authentication("TICKET_CREATE")

        """
        pass

class DoctestTicketProxy(object):
    r"""
    A proxy object for the ``ticket`` property of a
    :class:`DoctestServerProxy`.

    EXAMPLES::

        sage: from sage.dev.test.server_proxy import DoctestServerProxy
        sage: from sage.dev.test.trac_server import DoctestTracServer
        sage: DoctestServerProxy(DoctestTracServer()).ticket
        <sage.dev.test.server_proxy.DoctestTicketProxy at 0x...>

    """
    def __init__(self, server_proxy):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: type(trac._anonymous_server_proxy.ticket)
            <class 'sage.dev.test.server_proxy.DoctestTicketProxy'>

        """
        self._server_proxy = server_proxy

    def create(self, summary, description, attributes):
        r"""
        Create a new ticket and return its ticket number.

        EXAMPLES::

            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: trac._authenticated_server_proxy.ticket.create('summary', 'description', {})
            1

        """
        self._server_proxy._check_authentication("TICKET_CREATE")

        from trac_server import Ticket
        ticket = len(self._server_proxy._server.tickets)+1
        self._server_proxy._server.tickets[ticket] = Ticket(ticket, summary, description, attributes)
        return ticket

    def update(self, ticket, comment, attributes):
        r"""
        Add a ``comment`` and update ``attributes`` of ``ticket``.

        OUTPUT:

        Returns a fake URL of the ticket.

        EXAMPLES::

            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: ticket = trac._authenticated_server_proxy.ticket.create('summary', 'description', {})
            sage: trac._authenticated_server_proxy.ticket.update(ticket, 'comment', {'component':'algebra'})
            'https://trac.sagemath.org/ticket/1#comment:1'

        """
        self._server_proxy._check_authentication("TICKET_MODIFY")

        ticket = self._server_proxy._server.tickets[ticket]
        ticket.comments.append(comment)
        ticket.attributes = attributes

        from sage.env import TRAC_SERVER_URI
        import urlparse
        return urlparse.urljoin(TRAC_SERVER_URI, 'ticket/%s#comment:%s'%(ticket.id, len(ticket.comments)))

    def get(self, ticket):
        r"""
        Return a tuple ``(ticket, time_created, time_changed, attributes)`` for
        ``ticket``.

        EXAMPLES::

            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: ticket = trac._authenticated_server_proxy.ticket.create('summary', 'description', {})
            sage: trac._anonymous_server_proxy.ticket.get(ticket)
            [1, 'not implemented', 'not implemented', {'description': 'description', 'summary': 'summary'}]

        """
        ticket = self._server_proxy._server.tickets[ticket]
        return [ticket.id, ticket.time_created, ticket.time_changed, ticket.attributes]
