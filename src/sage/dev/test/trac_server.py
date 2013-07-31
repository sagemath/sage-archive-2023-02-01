r"""
Trac Server for Doctesting

This module provides a fake trac server which can be used for doctesting.

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

class DoctestTracServer(object):
    r"""
    A trac "server" which can be shared among instances of
    :class:`trac_interface.DoctestTracInterface` to keep track of tickets
    during doctesting.

    EXAMPLES::

        sage: from sage.dev.test.trac_server import DoctestTracServer
        sage: DoctestTracServer()
        <sage.dev.test.trac_server.DoctestTracServer object at 0x...>

    """
    def __init__(self):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: type(DoctestTracServer())
            <class 'sage.dev.test.trac_server.DoctestTracServer'>

        """
        self.tickets = {}

class Ticket(object):
    r"""
    Container for a ticket of a :class:`DoctestTracServer`.

    EXAMPLES::

        sage: from sage.dev.test.trac_server import Ticket
        sage: Ticket(1, "summary", "description", {})
        <sage.dev.test.trac_server.Ticket object at 0x...>

    """
    def __init__(self, ticket, summary, description, attributes):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.test.trac_server import Ticket
            sage: type(Ticket(1, "summary", "description", {}))
            <class 'sage.dev.test.trac_server.Ticket'>

        """
        attributes['summary'] = summary
        attributes['description'] = description
        self.id = ticket
        self.attributes = attributes
        self.time_created = 'not implemented'
        self.time_changed = 'not implemented'
        self.comments = []
