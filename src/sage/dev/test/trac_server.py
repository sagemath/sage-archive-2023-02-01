r"""
Trac Server for Doctesting

This module provides a fake trac server which can be used for doctesting.

AUTHORS:

- Julian Rueth: initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 Julian Rueth <julian.rueth@fsfe.org>
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
        from sage.dev.test.user_interface import DoctestUserInterface
        from sage.dev.test.config import DoctestConfig
        from sage.dev.sagedev import MASTER_BRANCH
        from sage.dev.git_interface import GitInterface
        config = DoctestConfig()
        self.git = GitInterface(config['git'], DoctestUserInterface(config['UI']))

        import os
        old_cwd = os.getcwd()
        os.chdir(config['git']['src'])
        try:
            self.git.super_silent.commit(allow_empty=True, message='initial commit')
            if MASTER_BRANCH != "master": self.git.super_silent.checkout("-b", MASTER_BRANCH)
            from sage.env import SAGE_VERSION
            self.git.super_silent.tag(SAGE_VERSION)
        finally:
            os.chdir(old_cwd)

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
        self.attachments = {}
