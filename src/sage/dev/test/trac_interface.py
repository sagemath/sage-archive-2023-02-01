r"""
Interface to trac Server for Doctesting

This module provides a subclass of
:class:`sage.dev.trac_interface.TracInterface` which can be used for
doctesting.

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
import sage.dev.trac_interface            


class DoctestTracInterface(sage.dev.trac_interface.TracInterface):
    r"""
    A :class:`sage.dev.trac_interface.TracInterface` which does not talk to the
    actual live server but to a :class:`sage.dev.test.trac_server.DoctestTracServer`.

    EXAMPLES::

        sage: from sage.dev.test.trac_interface import DoctestTracInterface
        sage: from sage.dev.test.config import DoctestConfig
        sage: from sage.dev.test.user_interface import DoctestUserInterface
        sage: from sage.dev.test.trac_server import DoctestTracServer
        sage: config = DoctestConfig()
        sage: UI = DoctestUserInterface(config['UI'])
        sage: DoctestTracInterface(config['trac'], UI, DoctestTracServer())
        <sage.dev.test.trac_interface.DoctestTracInterface object at 0x...>

    """
    def __init__(self, config, UI, server):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: UI = DoctestUserInterface(config['UI'])
            sage: type(DoctestTracInterface(config['trac'], UI, DoctestTracServer()))
            <class 'sage.dev.test.trac_interface.DoctestTracInterface'>
        """
        sage.dev.trac_interface.TracInterface.__init__(self, config, UI)
        self._server = server
        self._connected = True

    @property
    def _anonymous_server_proxy(self):
        r"""
        Return an non-authenticated proxy to the
        :class:`trac_server.DoctestTracServer` of this object.

        EXAMPLES::

            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: trac._anonymous_server_proxy
            <sage.dev.test.server_proxy.DoctestServerProxy object at 0x...>

        """
        if not self._connected:
            from sage.dev.trac_error import TracConnectionError
            raise TracConnectionError

        from server_proxy import DoctestServerProxy
        return DoctestServerProxy(self._server)

    @property
    def _authenticated_server_proxy(self):
        r"""
        Return an non-authenticated proxy to the
        :class:`trac_server.DoctestTracServer` of this object.

        EXAMPLES::

            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: trac._authenticated_server_proxy
            <sage.dev.test.server_proxy.AuthenticatedDoctestServerProxy object at 0x...>

        """
        if not self._connected:
            from sage.dev.trac_error import TracConnectionError
            raise TracConnectionError

        from server_proxy import AuthenticatedDoctestServerProxy
        return AuthenticatedDoctestServerProxy(self._server, self._username, self._password)
